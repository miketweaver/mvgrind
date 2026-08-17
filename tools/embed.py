#!/usr/bin/env python3
import sys

USAGE = "usage: ./embed.py <out.h> NAME=path [NAME=path ...]"


def literal(text: str) -> str:
    out = []
    for line in text.split("\n"):
        esc = line.replace("\\", "\\\\").replace('"', '\\"')
        out.append(f'    "{esc}\\n"')
    return "\n".join(out) if out else '    ""'


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(USAGE)
    out_path, specs = sys.argv[1], sys.argv[2:]

    chunks = [
        "#ifndef MV_KERNEL_SRC_H",
        "#define MV_KERNEL_SRC_H",
        "",
    ]
    for spec in specs:
        name, _, path = spec.partition("=")
        text = open(path, encoding="utf-8").read()
        chunks.append(f"static const char {name}[] =")
        chunks.append(literal(text) + ";")
        chunks.append("")
    chunks.append("#endif")

    open(out_path, "w", encoding="utf-8").write("\n".join(chunks) + "\n")
    print(f"wrote {out_path} ({len(specs)} sources embedded)")


if __name__ == "__main__":
    main()
