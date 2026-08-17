# mvgrind

GPU grinder for Meshtastic vanity NodeNums. Vendor-neutral OpenCL: one binary
runs on NVIDIA, AMD, Intel and POCL. A full 8-digit node id takes ~25 s on a
laptop RTX 3050 Ti (~172 M keys/s).

On PKI firmware builds, `nodenum = crc32(x25519_public_key)`. Both steps are
one-way, so a chosen id means searching the keyspace: pick a random clamped
scalar, step by 8 (one Edwards point addition per candidate instead of a full
ladder), check the CRC.

## Build

```sh
git clone --recursive https://github.com/miketweaver/mvgrind
cd mvgrind && make
make test
```

Needs a C compiler, `make`, `python3`, and your GPU vendor's OpenCL driver.
Already cloned without `--recursive`? Run `git submodule update --init --recursive`.

## Use

```sh
./mvgrind dc80                  # prefix: node id starting !dc80
./mvgrind dc801051              # exact id
./mvgrind 'dc80****'            # wildcard nibbles
./mvgrind 'dc80,801f,d0f0'      # a set: any match wins, no extra cost
./mvgrind --list                # enumerate OpenCL devices
./mvgrind --selftest            # CPU vectors + CPU/GPU differential
./mvgrind dc80 --bench 20       # measure throughput
./mvgrind -h                    # all options
```

Hits print to the terminal and append to `found.txt`. **That file is a private
key, guard it**. Want several ids? Grind them as one set; it's basically free.

The GPU self-test runs before every grind, and every hit is re-derived on the
CPU, so a broken kernel cannot produce a wrong key.

Verify a key independently (PyNaCl + zlib, shares no code with the grinder):

```sh
./tools/mvverify.py <private-key-hex-or-base64>
```

## Notes

- Keys are seeded from `getrandom(2)` per run: full 256-bit entropy, different
  key every run. Don't ever "simplify" this to a clock or counter seed; that
  makes every key brute-forceable.
- Only clamped keys are emitted. The firmware signs with a clamped copy of the
  scalar, so an unclamped key produces a node whose signatures don't verify.
- A 32-bit id is not an identity: anyone can grind a different key with the
  same id. It's a cosmetic label; security comes from the signature.

## License

GPL-3.0 (see `LICENSE`). No third-party source in-tree: the X25519 core is
extracted at build time from the pinned [Monocypher](https://github.com/LoupVaillant/Monocypher)
submodule (BSD-2-Clause OR CC0-1.0); [OpenCL-Headers](https://github.com/KhronosGroup/OpenCL-Headers)
(Apache-2.0) is used only if your distro has no `CL/cl.h`.
