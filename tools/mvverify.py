#!/usr/bin/env python3
import base64
import binascii
import sys
import zlib

from nacl.bindings import crypto_scalarmult_base

USAGE = """\
usage:
  ./mvverify.py <private-key-hex-or-base64>
  ./mvverify.py --selftest"""

CRC32_POLYNOMIAL = 0xEDB88320
CRC32_INITIAL = 0xFFFFFFFF


def crc32_erriez(buf: bytes) -> int:
    crc = CRC32_INITIAL
    for byte in buf:
        crc ^= byte
        for _ in range(8):
            mask = -(crc & 1) & 0xFFFFFFFF
            crc = (crc >> 1) ^ (CRC32_POLYNOMIAL & mask)
    return (~crc) & 0xFFFFFFFF


def is_clamped(sk: bytes) -> bool:
    return sk[0] & 0x07 == 0 and sk[31] & 0x80 == 0 and sk[31] & 0x40 == 0x40


def load_key(text: str) -> bytes:
    text = text.strip()
    try:
        raw = binascii.unhexlify(text)
        if len(raw) == 32:
            return raw
    except (binascii.Error, ValueError):
        pass
    raw = base64.b64decode(text)
    if len(raw) != 32:
        raise SystemExit(f"key must be 32 bytes, got {len(raw)}")
    return raw


def selftest() -> None:
    vectors = [b"", b"a", b"123456789", bytes(range(32)), b"\xff" * 32]
    for v in vectors:
        a, b = zlib.crc32(v) & 0xFFFFFFFF, crc32_erriez(v)
        assert a == b, f"CRC mismatch on {v!r}: zlib={a:08x} erriez={b:08x}"
    assert crc32_erriez(b"123456789") == 0xCBF43926, "check value wrong"
    print("selftest ok: zlib.crc32 == ErriezCRC32 transcription, check value 0xcbf43926")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(USAGE)
    if sys.argv[1] == "--selftest":
        selftest()
        return

    sk = load_key(sys.argv[1])
    pk = crypto_scalarmult_base(sk)

    zid = zlib.crc32(pk) & 0xFFFFFFFF
    eid = crc32_erriez(pk)
    if zid != eid:
        raise SystemExit(f"CRC implementations disagree: {zid:08x} vs {eid:08x}")

    print(f"private key : {sk.hex()}")
    print(f"public key  : {pk.hex()}")
    print(f"priv base64 : {base64.b64encode(sk).decode()}")
    print(f"pub  base64 : {base64.b64encode(pk).decode()}")
    print(f"node id     : !{zid:08x}")
    print(f"clamped     : {'yes' if is_clamped(sk) else 'NO, firmware XEdDSA will break'}")
    if not is_clamped(sk):
        sys.exit(2)


if __name__ == "__main__":
    main()
