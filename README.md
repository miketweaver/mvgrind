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

Needs a C compiler, `make`, `python3`, and your GPU vendor's OpenCL driver. On
macOS the system OpenCL framework is that driver, so nothing extra to install.
Already cloned without `--recursive`? Run `git submodule update --init --recursive`.

## Use

```sh
./mvgrind dc80                  # prefix: node id starting !dc80
./mvgrind dc801051              # exact id
./mvgrind 'dc80****'            # wildcard nibbles
./mvgrind 'dc80,801f,d0f0'      # a set: any match wins, no extra cost
./mvgrind --color crimson       # a node the apps paint crimson
./mvgrind --color '#dc143c'     # the same, by hex
./mvgrind dc --color crimson    # both: id starts !dc, and it's crimson
./mvgrind --color teal --tol 6  # near enough to teal, 2000x fewer keys
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

## Colors

Every client paints a node with the low 24 bits of its id, read straight as RGB.
[Android](https://github.com/meshtastic/Meshtastic-Android/blob/main/core/model/src/commonMain/kotlin/org/meshtastic/core/model/NodeColors.kt)
and [Apple](https://github.com/meshtastic/Meshtastic-Apple/blob/main/Meshtastic/Extensions/Color.swift)
agree on the formula. So `!8adc143c` is crimson, and picking a color is just
another pattern over the id, and `--color` compiles one for you.

`--color` takes `#rrggbb`, `#rgb`, or any of the 148 CSS color names. `--tol N`
widens it to anything within N per channel, which costs nothing extra to check
and finds a hit sooner:

| grind | keys (mean) | ~time at 172 M/s |
| --- | --- | --- |
| `--color crimson` | 16.8 M | 0.1 s |
| `--color crimson --tol 6` | 7.6 K | instant |
| `dc --color crimson` | 4.3 G | 25 s |
| `dc --color crimson --tol 6` | 2.0 M | 0.01 s |

The two patterns share bits: id nibbles 3-8 **are** the color channels, so
`dc80` already pins red to `0x80`. Ask for a color it can't be and mvgrind says
so up front rather than grinding forever:

```
$ ./mvgrind dc80 --color crimson
the id pattern and that color disagree on the red channel:
the pattern needs (byte & 0xff) == 0x80, the color needs 0xdc
```

Tolerance is capped by the 4096 targets the GPU holds: `--tol 7` is the widest
box that fits when all three channels are free, and much wider once an id
prefix pins one of them.

## Speed tricks

Grinding means making millions of keypairs a second and checking each one's id.
Three tricks get this from a crawl to ~172 M keys/s:

1. **Run it on the GPU.** A GPU has thousands of tiny cores. Each one grinds its
   own slice of the search at the same time, so we test thousands of keys in
   parallel instead of one at a time on the CPU.

2. **Don't build each key from scratch.** Normally every key needs a long, full
   calculation (~2550 steps of math). But we pick keys in order, and each key in
   the sequence differs from the last one by a fixed amount. So instead of
   rebuilding, we just take the previous key and nudge it forward, one cheap step
   (~9 steps of math) instead of 2550. This alone is roughly a 5x speedup.
   Source: the one-step "nudge" is
   [`mv_ge_add`](kernels/mv_edwards.h#L53) (a single elliptic-curve point
   addition).

3. **Share the most expensive step across a batch.** The one slow part left is a
   "division" that each key needs. There's a classic trick to do one division for
   a whole batch of keys instead of one per key. At a batch of 64, that expensive
   step nearly disappears, for another ~12x.
   Source: the batched-inversion loop (forward products, one `fe_invert`,
   backward pass) in
   [`mv_grind_batch`](kernels/mv_grind.cl#L82-L105).

Together: ~690x faster than the CPU. `--ladder` runs the old slow way and
`--batch 1` turns off trick 3, so you can measure each step yourself with
`--bench`.

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
