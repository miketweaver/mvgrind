#!/usr/bin/env python3
import re
import sys

USAGE = "usage: ./extract_monocypher.py <monocypher.c> <out.inc>"

START_BANNER = "/// Arithmetic modulo 2^255 - 19 ///"
END_BANNER = "/// Arithmetic modulo L ///"

SQRTM1_LOCAL = """\tconst fe sqrtm1 = {
\t\t-32595792, -7943725, 9377950, 3500415, 12389472,
\t\t-272473, -25146209, -2005654, 326686, 11406482,
\t};
"""

HEADER = """\
/* X25519 core, extracted from Monocypher by tools/extract_monocypher.py.
 * DO NOT EDIT BY HAND, rerun the extractor instead.
 *
 * Upstream: https://github.com/LoupVaillant/Monocypher
 * Copyright (c) 2017-2020, Loup Vaillant
 * SPDX-License-Identifier: BSD-2-Clause OR CC0-1.0
 */

"""


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(USAGE)
    src = open(sys.argv[1], encoding="utf-8").read()

    start = src.index(START_BANNER)
    start = src.index("\n", src.index("////", start + len(START_BANNER))) + 1
    end = src.rindex("///////", 0, src.index(END_BANNER))
    body = src[start:end]

    body, n = re.subn(
        r"^static const fe [A-Za-z0-9_]+\s*=\s*\{.*?\};\n",
        "",
        body,
        flags=re.DOTALL | re.MULTILINE,
    )
    if n != 8:
        raise SystemExit(f"expected 8 constant tables, removed {n}, upstream changed")

    marker = "static int invsqrt(fe isr, const fe x)\n{\n"
    if marker not in body:
        raise SystemExit("invsqrt signature changed, update the extractor")
    body = body.replace(marker, marker + SQRTM1_LOCAL, 1)

    body, n_bp = re.subn(
        r"\bstatic const u8 base_point\[32\] = \{9\};",
        "const u8 base_point[32] = {9};",
        body,
    )
    if n_bp != 1:
        raise SystemExit(f"expected 1 base_point static, rewrote {n_bp}")

    n_static = 0
    for fn in ("crypto_eddsa_trim_scalar", "crypto_x25519_public_key", "crypto_x25519"):
        body, k = re.subn(rf"^void {fn}\(", f"static void {fn}(", body, flags=re.MULTILINE)
        n_static += k
    if n_static != 3:
        raise SystemExit(f"expected 3 exported functions to hide, got {n_static}")

    if "static const fe" in body:
        raise SystemExit("a program-scope fe constant survived, not OpenCL-legal")
    leftover = re.search(r"^[ \t]+static\b", body, flags=re.MULTILINE)
    if leftover:
        raise SystemExit(f"function-scope static survived: {leftover.group(0)!r}")
    out = HEADER + body.strip() + "\n"
    open(sys.argv[2], "w", encoding="utf-8").write(out)
    print(f"wrote {sys.argv[2]}: {out.count(chr(10))} lines, {n} constant tables removed")


if __name__ == "__main__":
    main()
