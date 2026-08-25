#define CL_TARGET_OPENCL_VERSION 120
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl.h>
#include <errno.h>
#include <inttypes.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#if defined(__has_include)
#if __has_include(<sys/random.h>)
#include <sys/random.h>
#define MV_HAVE_GETRANDOM 1
#endif
#endif

#include "kernels/mv_prelude.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "monocypher_x25519.inc"
#pragma GCC diagnostic pop
#include "kernels/mv_edwards.h"
#include "mv_kernel_src.h"

#define MV_GE_INTS 40
#define MV_MAX_DEVICES 16
#define MV_OUT_CAP 256
#define MV_MAX_TARGETS 4096

/* ------------------------------------------------------------------------- */
/* Incremental (Edwards) walk: secrets are base + 8k, so consecutive public   */
/* keys differ by the constant point 8G.                                      */
/* ------------------------------------------------------------------------- */

static void mv_secret_at(uint8_t out[32], const uint8_t base[32], uint64_t k)
{
    COPY(out, base, 32);
    mv_add8k(out, k);
}

static void scalar_from_u64(uint8_t out[32], uint64_t v)
{
    ZERO(out, 32);
    for (int i = 0; i < 8; i++)
        out[i] = (uint8_t)(v >> (8 * i));
}

static int mv_build_walk(int32_t *p_start, int32_t *jumps, unsigned n_jumps, int32_t *step, const uint8_t base[32],
                         uint64_t k_base, uint32_t iters)
{
    if (!p_start || !step || (n_jumps && !jumps) || iters == 0)
        return -1;

    uint8_t s[32];
    mv_ge p;
    mv_secret_at(s, base, k_base);
    mv_ge_scalarmult_base(&p, s);
    memcpy(p_start, &p, sizeof(mv_ge));

    mv_ge g, eight_g;
    mv_ge_cached c;
    mv_ge_base(&g);
    mv_ge_double(&eight_g, &g);
    mv_ge_double(&eight_g, &eight_g);
    mv_ge_double(&eight_g, &eight_g);
    mv_ge_cache(&c, &eight_g);
    memcpy(step, &c, sizeof(mv_ge_cached));

    uint8_t js[32];
    mv_ge jp;
    scalar_from_u64(js, (uint64_t)iters * 8u);
    mv_ge_scalarmult_base(&jp, js);
    for (unsigned j = 0; j < n_jumps; j++) {
        mv_ge_cache(&c, &jp);
        memcpy(jumps + (size_t)j * MV_GE_INTS, &c, sizeof(mv_ge_cached));
        mv_ge_double(&jp, &jp);
    }

    WIPE_BUFFER(s);
    return 0;
}

/* ------------------------------------------------------------------------- */
/* Patterns: a mask plus a set of accepted masked values                      */
/* ------------------------------------------------------------------------- */

typedef struct {
    uint32_t mask;
    uint32_t *targets;
    uint32_t n_targets;
} mv_pattern;

static int parse_one(const char *spec, uint32_t *val, uint32_t *mask, char *msg, size_t msglen)
{
    if (*spec == '!')
        spec++;
    size_t n = strlen(spec);
    if (n == 0 || n > 8) {
        if (msg)
            snprintf(msg, msglen, "'%s': need 1-8 hex digits or wildcards", spec);
        return -1;
    }
    uint32_t v = 0, m = 0;
    for (size_t i = 0; i < n; i++) {
        char c = spec[i];
        v <<= 4;
        m <<= 4;
        if (c == '*' || c == '?' || c == '.' || c == 'x' || c == 'X')
            continue;
        int d;
        if (c >= '0' && c <= '9')
            d = c - '0';
        else if (c >= 'a' && c <= 'f')
            d = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F')
            d = c - 'A' + 10;
        else {
            if (msg)
                snprintf(msg, msglen, "'%s': bad character '%c'", spec, c);
            return -1;
        }
        v |= (uint32_t)d;
        m |= 0xFu;
    }
    int shift = (int)(8 - n) * 4;
    *val = v << shift;
    *mask = m << shift;
    return 0;
}

static int mv_pattern_parse(mv_pattern *p, const char *spec, char *msg, size_t msglen)
{
    memset(p, 0, sizeof(*p));
    if (!spec || !*spec) {
        if (msg)
            snprintf(msg, msglen, "empty pattern");
        return -1;
    }

    uint32_t count = 1;
    for (const char *q = spec; *q; q++)
        if (*q == ',')
            count++;

    p->targets = (uint32_t *)calloc(count, sizeof(uint32_t));
    char *dup = strdup(spec);
    if (!p->targets || !dup) {
        free(p->targets);
        free(dup);
        p->targets = NULL;
        if (msg)
            snprintf(msg, msglen, "out of memory");
        return -1;
    }

    int ok = 0;
    char *save = NULL;
    for (char *tok = strtok_r(dup, ",", &save); tok; tok = strtok_r(NULL, ",", &save)) {
        while (*tok == ' ')
            tok++;
        uint32_t v, m;
        if (parse_one(tok, &v, &m, msg, msglen) != 0)
            goto fail;
        if (p->n_targets == 0) {
            p->mask = m;
        } else if (m != p->mask) {
            if (msg)
                snprintf(msg, msglen,
                         "all patterns in a set must fix the same nibbles "
                         "(got mask %08x then %08x)",
                         p->mask, m);
            goto fail;
        }
        p->targets[p->n_targets++] = v & m;
    }
    if (p->n_targets == 0) {
        if (msg)
            snprintf(msg, msglen, "no usable pattern terms");
        goto fail;
    }
    ok = 1;

fail:
    free(dup);
    if (!ok) {
        free(p->targets);
        p->targets = NULL;
        return -1;
    }
    return 0;
}

static int mv_pattern_match(const mv_pattern *p, uint32_t nodenum)
{
    uint32_t v = nodenum & p->mask;
    for (uint32_t i = 0; i < p->n_targets; i++)
        if (p->targets[i] == v)
            return 1;
    return 0;
}

static void mv_pattern_free(mv_pattern *p)
{
    free(p->targets);
    p->targets = NULL;
    p->n_targets = 0;
}

static double mv_pattern_expected_keys(const mv_pattern *p)
{
    int bits = 0;
    for (uint32_t m = p->mask; m; m >>= 1)
        bits += (int)(m & 1u);
    double space = 1.0;
    for (int i = 0; i < bits; i++)
        space *= 2.0;
    return space / (p->n_targets ? (double)p->n_targets : 1.0);
}

/* ------------------------------------------------------------------------- */
/* Colors: every client paints a node with the low 24 bits of its NodeNum as  */
/* RGB, so a chosen color is just a pattern over those bits.                  */
/* ------------------------------------------------------------------------- */

#define MV_COLOR_MASK 0x00FFFFFFu

typedef struct {
    uint8_t r, g, b;
} mv_color;

/* The 148 CSS Color Module Level 4 keywords, sorted for bsearch. */
typedef struct {
    const char *name;
    uint32_t rgb;
} mv_css_color;

static const mv_css_color mv_css_colors[] = {
    {"aliceblue", 0xf0f8ff},
    {"antiquewhite", 0xfaebd7},
    {"aqua", 0x00ffff},
    {"aquamarine", 0x7fffd4},
    {"azure", 0xf0ffff},
    {"beige", 0xf5f5dc},
    {"bisque", 0xffe4c4},
    {"black", 0x000000},
    {"blanchedalmond", 0xffebcd},
    {"blue", 0x0000ff},
    {"blueviolet", 0x8a2be2},
    {"brown", 0xa52a2a},
    {"burlywood", 0xdeb887},
    {"cadetblue", 0x5f9ea0},
    {"chartreuse", 0x7fff00},
    {"chocolate", 0xd2691e},
    {"coral", 0xff7f50},
    {"cornflowerblue", 0x6495ed},
    {"cornsilk", 0xfff8dc},
    {"crimson", 0xdc143c},
    {"cyan", 0x00ffff},
    {"darkblue", 0x00008b},
    {"darkcyan", 0x008b8b},
    {"darkgoldenrod", 0xb8860b},
    {"darkgray", 0xa9a9a9},
    {"darkgreen", 0x006400},
    {"darkgrey", 0xa9a9a9},
    {"darkkhaki", 0xbdb76b},
    {"darkmagenta", 0x8b008b},
    {"darkolivegreen", 0x556b2f},
    {"darkorange", 0xff8c00},
    {"darkorchid", 0x9932cc},
    {"darkred", 0x8b0000},
    {"darksalmon", 0xe9967a},
    {"darkseagreen", 0x8fbc8f},
    {"darkslateblue", 0x483d8b},
    {"darkslategray", 0x2f4f4f},
    {"darkslategrey", 0x2f4f4f},
    {"darkturquoise", 0x00ced1},
    {"darkviolet", 0x9400d3},
    {"deeppink", 0xff1493},
    {"deepskyblue", 0x00bfff},
    {"dimgray", 0x696969},
    {"dimgrey", 0x696969},
    {"dodgerblue", 0x1e90ff},
    {"firebrick", 0xb22222},
    {"floralwhite", 0xfffaf0},
    {"forestgreen", 0x228b22},
    {"fuchsia", 0xff00ff},
    {"gainsboro", 0xdcdcdc},
    {"ghostwhite", 0xf8f8ff},
    {"gold", 0xffd700},
    {"goldenrod", 0xdaa520},
    {"gray", 0x808080},
    {"green", 0x008000},
    {"greenyellow", 0xadff2f},
    {"grey", 0x808080},
    {"honeydew", 0xf0fff0},
    {"hotpink", 0xff69b4},
    {"indianred", 0xcd5c5c},
    {"indigo", 0x4b0082},
    {"ivory", 0xfffff0},
    {"khaki", 0xf0e68c},
    {"lavender", 0xe6e6fa},
    {"lavenderblush", 0xfff0f5},
    {"lawngreen", 0x7cfc00},
    {"lemonchiffon", 0xfffacd},
    {"lightblue", 0xadd8e6},
    {"lightcoral", 0xf08080},
    {"lightcyan", 0xe0ffff},
    {"lightgoldenrodyellow", 0xfafad2},
    {"lightgray", 0xd3d3d3},
    {"lightgreen", 0x90ee90},
    {"lightgrey", 0xd3d3d3},
    {"lightpink", 0xffb6c1},
    {"lightsalmon", 0xffa07a},
    {"lightseagreen", 0x20b2aa},
    {"lightskyblue", 0x87cefa},
    {"lightslategray", 0x778899},
    {"lightslategrey", 0x778899},
    {"lightsteelblue", 0xb0c4de},
    {"lightyellow", 0xffffe0},
    {"lime", 0x00ff00},
    {"limegreen", 0x32cd32},
    {"linen", 0xfaf0e6},
    {"magenta", 0xff00ff},
    {"maroon", 0x800000},
    {"mediumaquamarine", 0x66cdaa},
    {"mediumblue", 0x0000cd},
    {"mediumorchid", 0xba55d3},
    {"mediumpurple", 0x9370db},
    {"mediumseagreen", 0x3cb371},
    {"mediumslateblue", 0x7b68ee},
    {"mediumspringgreen", 0x00fa9a},
    {"mediumturquoise", 0x48d1cc},
    {"mediumvioletred", 0xc71585},
    {"midnightblue", 0x191970},
    {"mintcream", 0xf5fffa},
    {"mistyrose", 0xffe4e1},
    {"moccasin", 0xffe4b5},
    {"navajowhite", 0xffdead},
    {"navy", 0x000080},
    {"oldlace", 0xfdf5e6},
    {"olive", 0x808000},
    {"olivedrab", 0x6b8e23},
    {"orange", 0xffa500},
    {"orangered", 0xff4500},
    {"orchid", 0xda70d6},
    {"palegoldenrod", 0xeee8aa},
    {"palegreen", 0x98fb98},
    {"paleturquoise", 0xafeeee},
    {"palevioletred", 0xdb7093},
    {"papayawhip", 0xffefd5},
    {"peachpuff", 0xffdab9},
    {"peru", 0xcd853f},
    {"pink", 0xffc0cb},
    {"plum", 0xdda0dd},
    {"powderblue", 0xb0e0e6},
    {"purple", 0x800080},
    {"rebeccapurple", 0x663399},
    {"red", 0xff0000},
    {"rosybrown", 0xbc8f8f},
    {"royalblue", 0x4169e1},
    {"saddlebrown", 0x8b4513},
    {"salmon", 0xfa8072},
    {"sandybrown", 0xf4a460},
    {"seagreen", 0x2e8b57},
    {"seashell", 0xfff5ee},
    {"sienna", 0xa0522d},
    {"silver", 0xc0c0c0},
    {"skyblue", 0x87ceeb},
    {"slateblue", 0x6a5acd},
    {"slategray", 0x708090},
    {"slategrey", 0x708090},
    {"snow", 0xfffafa},
    {"springgreen", 0x00ff7f},
    {"steelblue", 0x4682b4},
    {"tan", 0xd2b48c},
    {"teal", 0x008080},
    {"thistle", 0xd8bfd8},
    {"tomato", 0xff6347},
    {"turquoise", 0x40e0d0},
    {"violet", 0xee82ee},
    {"wheat", 0xf5deb3},
    {"white", 0xffffff},
    {"whitesmoke", 0xf5f5f5},
    {"yellow", 0xffff00},
    {"yellowgreen", 0x9acd32},
};

static int mv_hexdigit(char c)
{
    if (c >= '0' && c <= '9')
        return c - '0';
    if (c >= 'a' && c <= 'f')
        return c - 'a' + 10;
    if (c >= 'A' && c <= 'F')
        return c - 'A' + 10;
    return -1;
}

/* Fold to the table's spelling: lowercase, letters and digits only, so
 * "Dark Slate Blue", "dark-slate-blue" and "darkslateblue" all land. */
static int mv_color_fold(char *dst, size_t cap, const char *src)
{
    size_t o = 0;
    for (; *src; src++) {
        char c = *src;
        if (c == ' ' || c == '-' || c == '_')
            continue;
        if (c >= 'A' && c <= 'Z')
            c = (char)(c - 'A' + 'a');
        if (!((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9')))
            return -1;
        if (o + 1 >= cap)
            return -1;
        dst[o++] = c;
    }
    dst[o] = 0;
    return o ? 0 : -1;
}

static int mv_color_cmp(const void *key, const void *ent)
{
    return strcmp((const char *)key, ((const mv_css_color *)ent)->name);
}

static mv_color mv_color_from_u32(uint32_t rgb)
{
    mv_color c;
    c.r = (uint8_t)(rgb >> 16);
    c.g = (uint8_t)(rgb >> 8);
    c.b = (uint8_t)rgb;
    return c;
}

static uint32_t mv_color_to_u32(mv_color c)
{
    return ((uint32_t)c.r << 16) | ((uint32_t)c.g << 8) | c.b;
}

/* Accepts "#rrggbb", "#rgb", the bare hex forms, and any CSS color keyword. */
static int mv_color_parse(const char *spec, mv_color *out, char *msg, size_t msglen)
{
    if (!spec || !*spec) {
        if (msg)
            snprintf(msg, msglen, "empty color");
        return -1;
    }
    while (*spec == ' ')
        spec++;

    int forced_hex = (*spec == '#');
    const char *h = forced_hex ? spec + 1 : spec;
    size_t n = strlen(h);

    if (!forced_hex) {
        char fold[64];
        if (mv_color_fold(fold, sizeof fold, spec) == 0) {
            const void *hit = bsearch(fold, mv_css_colors, sizeof(mv_css_colors) / sizeof(mv_css_colors[0]),
                                      sizeof(mv_css_colors[0]), mv_color_cmp);
            if (hit) {
                *out = mv_color_from_u32(((const mv_css_color *)hit)->rgb);
                return 0;
            }
        }
    }

    if (n != 3 && n != 6) {
        if (msg)
            snprintf(msg, msglen, "'%s': not a CSS color name, and hex needs 3 or 6 digits", spec);
        return -1;
    }
    int d[6];
    for (size_t i = 0; i < n; i++) {
        d[i] = mv_hexdigit(h[i]);
        if (d[i] < 0) {
            if (msg)
                snprintf(msg, msglen, "'%s': bad hex digit '%c'", spec, h[i]);
            return -1;
        }
    }
    if (n == 3) {
        out->r = (uint8_t)(d[0] * 17);
        out->g = (uint8_t)(d[1] * 17);
        out->b = (uint8_t)(d[2] * 17);
    } else {
        out->r = (uint8_t)(d[0] * 16 + d[1]);
        out->g = (uint8_t)(d[2] * 16 + d[3]);
        out->b = (uint8_t)(d[4] * 16 + d[5]);
    }
    return 0;
}

static mv_color mv_color_of_nodenum(uint32_t nodenum)
{
    return mv_color_from_u32(nodenum & MV_COLOR_MASK);
}

/* "#dc143c", optionally followed by a truecolor swatch when stdout is a tty. */
static void mv_color_fmt(char *dst, size_t cap, mv_color c, int swatch)
{
    if (swatch)
        snprintf(dst, cap, "#%02x%02x%02x  \033[48;2;%u;%u;%um    \033[0m", c.r, c.g, c.b, c.r, c.g, c.b);
    else
        snprintf(dst, cap, "#%02x%02x%02x", c.r, c.g, c.b);
}

static const char *const mv_chan_name[3] = {"red", "green", "blue"};
static const int mv_chan_shift[3] = {16, 8, 0};

/* Values in lo..hi that also satisfy the id pattern's bits for this channel. */
static int mv_chan_allowed(uint8_t *out, int lo, int hi, uint32_t pmask, uint32_t pval)
{
    int n = 0;
    for (int x = lo; x <= hi; x++)
        if (((uint32_t)x & pmask) == (pval & pmask))
            out[n++] = (uint8_t)x;
    return n;
}

/* Fold a color box (every RGB within +/-tol of c, clamped to 0..255) into
 * *p, which may already hold an id pattern (and may be empty). The two can
 * overlap: id nibbles 3-8 are the color channels, so 'dc80' pins red to 0x80. */
static int mv_pattern_apply_color(mv_pattern *p, mv_color c, int tol, char *msg, size_t msglen)
{
    if (tol < 0 || tol > 255) {
        if (msg)
            snprintf(msg, msglen, "tolerance must be 0-255 (got %d)", tol);
        return -1;
    }

    const uint8_t want[3] = {c.r, c.g, c.b};
    int lo[3], hi[3];
    for (int k = 0; k < 3; k++) {
        lo[k] = want[k] - tol < 0 ? 0 : want[k] - tol;
        hi[k] = want[k] + tol > 255 ? 255 : want[k] + tol;
    }

    const uint32_t id_mask = p->n_targets ? p->mask : 0;
    const uint32_t none = 0;
    const uint32_t *ids = p->n_targets ? p->targets : &none;
    const uint32_t n_ids = p->n_targets ? p->n_targets : 1;

    uint32_t *out = (uint32_t *)calloc(MV_MAX_TARGETS, sizeof(uint32_t));
    if (!out) {
        if (msg)
            snprintf(msg, msglen, "out of memory");
        return -1;
    }

    uint32_t n_out = 0;
    int chan_ok[3] = {0, 0, 0};
    for (uint32_t i = 0; i < n_ids; i++) {
        uint8_t a[3][256];
        int na[3];
        uint64_t prod = 1;
        for (int k = 0; k < 3; k++) {
            const uint32_t pm = (id_mask >> mv_chan_shift[k]) & 0xFFu;
            const uint32_t pv = (ids[i] >> mv_chan_shift[k]) & 0xFFu;
            na[k] = mv_chan_allowed(a[k], lo[k], hi[k], pm, pv);
            if (na[k])
                chan_ok[k] = 1;
            prod *= (uint64_t)na[k];
        }
        if (prod == 0)
            continue;
        if (n_out + prod > MV_MAX_TARGETS) {
            if (msg)
                snprintf(msg, msglen,
                         "--tol %d needs more than %d targets, which is all the GPU "
                         "will hold; use a tighter tolerance",
                         tol, MV_MAX_TARGETS);
            free(out);
            return -1;
        }
        const uint32_t keep = ids[i] & ~MV_COLOR_MASK;
        for (int ri = 0; ri < na[0]; ri++)
            for (int gi = 0; gi < na[1]; gi++)
                for (int bi = 0; bi < na[2]; bi++)
                    out[n_out++] = keep | ((uint32_t)a[0][ri] << 16) | ((uint32_t)a[1][gi] << 8) | a[2][bi];
    }

    if (n_out == 0) {
        int k = 0;
        while (k < 3 && chan_ok[k])
            k++;
        if (msg && k < 3) {
            const uint32_t pm = (id_mask >> mv_chan_shift[k]) & 0xFFu;
            char wanted[16];
            if (lo[k] == hi[k])
                snprintf(wanted, sizeof wanted, "0x%02x", (unsigned)lo[k]);
            else
                snprintf(wanted, sizeof wanted, "0x%02x-0x%02x", (unsigned)lo[k], (unsigned)hi[k]);
            snprintf(msg, msglen,
                     "the id pattern and that color disagree on the %s channel: "
                     "the pattern needs (byte & 0x%02x) == 0x%02x, the color needs %s",
                     mv_chan_name[k], pm, (ids[0] >> mv_chan_shift[k]) & pm, wanted);
        } else if (msg) {
            snprintf(msg, msglen, "no node id satisfies both the id pattern and that color");
        }
        free(out);
        return -1;
    }

    free(p->targets);
    p->targets = out;
    p->n_targets = n_out;
    p->mask = id_mask | MV_COLOR_MASK;
    return 0;
}

#define MV_CFAIL(...)                                                                                                            \
    do {                                                                                                                         \
        snprintf(err, errlen, __VA_ARGS__);                                                                                      \
        mv_pattern_free(&p);                                                                                                     \
        return -1;                                                                                                               \
    } while (0)

static int mv_selftest_color(char *err, size_t errlen)
{
    mv_pattern p;
    memset(&p, 0, sizeof p);
    char msg[256];

    static const struct {
        const char *spec;
        uint32_t rgb;
    } good[] = {
        {"crimson", 0xdc143c},       {"CRIMSON", 0xdc143c}, {"Dark Slate Blue", 0x483d8b}, {"dark-slate_blue", 0x483d8b},
        {"rebeccapurple", 0x663399}, {"#dc143c", 0xdc143c}, {"dc143c", 0xdc143c},          {"#f00", 0xff0000},
        {"0f8", 0x00ff88},           {"#000", 0x000000},    {"white", 0xffffff},
    };
    for (size_t i = 0; i < sizeof(good) / sizeof(good[0]); i++) {
        mv_color c;
        if (mv_color_parse(good[i].spec, &c, msg, sizeof msg) != 0)
            MV_CFAIL("color '%s' should parse (%s)", good[i].spec, msg);
        if (mv_color_to_u32(c) != good[i].rgb)
            MV_CFAIL("color '%s' gave %06x, want %06x", good[i].spec, mv_color_to_u32(c), good[i].rgb);
    }

    static const char *const bad[] = {"", "notacolor", "#12345", "#gg0000", "#dc143c1", "12 34"};
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
        mv_color c;
        if (mv_color_parse(bad[i], &c, msg, sizeof msg) == 0)
            MV_CFAIL("color '%s' should have been rejected", bad[i]);
    }

    if (mv_color_to_u32(mv_color_of_nodenum(0xabdc143cu)) != 0xdc143c)
        MV_CFAIL("nodenum -> color is not the low 24 bits");

    /* Color alone: the low 24 bits, nothing else. */
    const mv_color crimson = mv_color_from_u32(0xdc143c);
    if (mv_pattern_apply_color(&p, crimson, 0, msg, sizeof msg) != 0)
        MV_CFAIL("exact color should compile (%s)", msg);
    if (p.mask != MV_COLOR_MASK || p.n_targets != 1 || p.targets[0] != 0x00dc143cu)
        MV_CFAIL("exact color gave mask %08x / %u targets", p.mask, p.n_targets);
    if (!mv_pattern_match(&p, 0xffdc143cu) || mv_pattern_match(&p, 0xffdc143du))
        MV_CFAIL("exact color matches the wrong ids");
    mv_pattern_free(&p);

    /* A tolerance box is (2t+1)^3 targets, clamped at the ends of each channel. */
    if (mv_pattern_apply_color(&p, crimson, 1, msg, sizeof msg) != 0 || p.n_targets != 27)
        MV_CFAIL("+/-1 box should be 27 targets, got %u", p.n_targets);
    mv_pattern_free(&p);
    if (mv_pattern_apply_color(&p, mv_color_from_u32(0x000000), 2, msg, sizeof msg) != 0 || p.n_targets != 27)
        MV_CFAIL("black +/-2 should clamp to 27 targets, got %u", p.n_targets);
    mv_pattern_free(&p);
    if (mv_pattern_apply_color(&p, mv_color_from_u32(0xffffff), 2, msg, sizeof msg) != 0 || p.n_targets != 27)
        MV_CFAIL("white +/-2 should clamp to 27 targets, got %u", p.n_targets);
    mv_pattern_free(&p);
    if (mv_pattern_apply_color(&p, crimson, 8, msg, sizeof msg) == 0)
        MV_CFAIL("+/-8 (4913 targets) should exceed the target cap");
    mv_pattern_free(&p);

    /* An id prefix short enough to miss the color channels combines cleanly. */
    if (mv_pattern_parse(&p, "dc", msg, sizeof msg) != 0)
        MV_CFAIL("pattern 'dc' should parse (%s)", msg);
    if (mv_pattern_apply_color(&p, crimson, 0, msg, sizeof msg) != 0)
        MV_CFAIL("'dc' + crimson should compile (%s)", msg);
    if (p.mask != 0xFFFFFFFFu || p.n_targets != 1 || p.targets[0] != 0xdcdc143cu)
        MV_CFAIL("'dc' + crimson gave mask %08x / %u targets / %08x", p.mask, p.n_targets, p.targets[0]);
    mv_pattern_free(&p);

    /* Nibbles 3-8 are the color channels, so a longer prefix can contradict it. */
    if (mv_pattern_parse(&p, "dc80", msg, sizeof msg) != 0)
        MV_CFAIL("pattern 'dc80' should parse (%s)", msg);
    if (mv_pattern_apply_color(&p, crimson, 0, msg, sizeof msg) == 0)
        MV_CFAIL("'dc80' pins red to 0x80 and must not accept crimson");
    mv_pattern_free(&p);

    /* ...and agree when it matches. */
    if (mv_pattern_parse(&p, "dcdc", msg, sizeof msg) != 0)
        MV_CFAIL("pattern 'dcdc' should parse (%s)", msg);
    if (mv_pattern_apply_color(&p, crimson, 0, msg, sizeof msg) != 0)
        MV_CFAIL("'dcdc' + crimson should compile (%s)", msg);
    if (p.n_targets != 1 || p.targets[0] != 0xdcdc143cu)
        MV_CFAIL("'dcdc' + crimson gave %u targets / %08x", p.n_targets, p.targets[0]);
    mv_pattern_free(&p);

    /* Tolerance survives a pinned channel: red fixed, green and blue free. */
    if (mv_pattern_parse(&p, "dcdc", msg, sizeof msg) != 0)
        MV_CFAIL("pattern 'dcdc' should parse (%s)", msg);
    if (mv_pattern_apply_color(&p, crimson, 8, msg, sizeof msg) != 0)
        MV_CFAIL("'dcdc' + crimson +/-8 should compile (%s)", msg);
    if (p.n_targets != 17 * 17)
        MV_CFAIL("'dcdc' + crimson +/-8 should be 289 targets, got %u", p.n_targets);
    mv_pattern_free(&p);

    return 0;
}

#undef MV_CFAIL

/* ------------------------------------------------------------------------- */
/* CPU self-test                                                              */
/* ------------------------------------------------------------------------- */

static int hexbytes(uint8_t *out, size_t n, const char *hex)
{
    for (size_t i = 0; i < n; i++) {
        unsigned v;
        if (sscanf(hex + 2 * i, "%2x", &v) != 1)
            return -1;
        out[i] = (uint8_t)v;
    }
    return 0;
}

static int mv_selftest_cpu(char *err, size_t errlen)
{
    struct {
        const char *sk, *pk;
    } v[] = {
        {"77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a",
         "8520f0098930a754748b7ddcb43ef75a0dbf3a0d26381af4eba4a98eaa9b4e6a"},
        {"5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb",
         "de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f"},
    };
    for (size_t i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
        uint8_t sk[32], want[32], got[32];
        if (hexbytes(sk, 32, v[i].sk) || hexbytes(want, 32, v[i].pk)) {
            snprintf(err, errlen, "bad builtin test vector");
            return -1;
        }
        crypto_x25519_public_key(got, sk);
        if (memcmp(got, want, 32) != 0) {
            snprintf(err, errlen, "RFC 7748 vector %zu failed", i);
            return -1;
        }
    }

    if (mv_crc32_core((const u8 *)"123456789", 9) != 0xCBF43926u) {
        snprintf(err, errlen, "CRC-32 check value wrong (expected 0xcbf43926)");
        return -1;
    }

    uint8_t s[32];
    memset(s, 0xFF, 32);
    mv_clamp(s);
    if (!mv_is_clamped(s)) {
        snprintf(err, errlen, "clamp/is_clamped disagree");
        return -1;
    }
    return mv_selftest_color(err, errlen);
}

/* ------------------------------------------------------------------------- */
/* OpenCL host                                                                */
/* ------------------------------------------------------------------------- */

typedef struct {
    uint8_t secret[32];
    uint8_t pub[32];
    uint32_t nodenum;
    int device_index;
} mv_hit;

typedef struct {
    char platform[128];
    char name[128];
    unsigned compute_units;
    unsigned clock_mhz;
    unsigned long long global_mem;
    int is_gpu;
} mv_device_info;

typedef struct {
    int device_index;
    unsigned iters_per_thread;
    size_t global_size;
    const char *kernel_dir;
    int use_ladder;
    unsigned batch;
} mv_grind_opts;

typedef int (*mv_hit_cb)(const mv_hit *hit, void *user);
typedef void (*mv_progress_cb)(uint64_t keys_tried, double keys_per_sec, void *user);

struct mv_dev {
    cl_device_id id;
    cl_context ctx;
    cl_command_queue q;
    cl_program prog;
    cl_kernel k_grind;
    cl_kernel k_batch;
    cl_kernel k_derive;
    cl_mem m_base, m_targets, m_count, m_k, m_id;
    cl_mem m_pstart, m_jumps, m_step;
    size_t global_size;
    size_t ramp;
    unsigned iters;
    unsigned n_jumps;
    int index;
    char name[128];
};

typedef struct {
    struct mv_dev dev[MV_MAX_DEVICES];
    int n_dev;
    mv_grind_opts opts;
} mv_ctx;

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

#define CHECK(e, what)                                                                                                           \
    do {                                                                                                                         \
        cl_int _e = (e);                                                                                                         \
        if (_e != CL_SUCCESS) {                                                                                                  \
            snprintf(err, errlen, "%s: OpenCL error %d", (what), (int)_e);                                                       \
            return -1;                                                                                                           \
        }                                                                                                                        \
    } while (0)

struct raw_dev {
    cl_platform_id plat;
    cl_device_id dev;
};

static int collect_devices(struct raw_dev *out, int cap)
{
    cl_uint nplat = 0;
    if (clGetPlatformIDs(0, NULL, &nplat) != CL_SUCCESS || nplat == 0)
        return 0;
    cl_platform_id plats[16];
    if (nplat > 16)
        nplat = 16;
    clGetPlatformIDs(nplat, plats, NULL);

    int n = 0;
    for (cl_uint i = 0; i < nplat && n < cap; i++) {
        cl_uint nd = 0;
        if (clGetDeviceIDs(plats[i], CL_DEVICE_TYPE_ALL, 0, NULL, &nd) != CL_SUCCESS || !nd)
            continue;
        cl_device_id ds[MV_MAX_DEVICES];
        if (nd > MV_MAX_DEVICES)
            nd = MV_MAX_DEVICES;
        clGetDeviceIDs(plats[i], CL_DEVICE_TYPE_ALL, nd, ds, NULL);
        for (cl_uint j = 0; j < nd && n < cap; j++) {
            out[n].plat = plats[i];
            out[n].dev = ds[j];
            n++;
        }
    }
    return n;
}

static int mv_list_devices(mv_device_info *buf, int cap)
{
    struct raw_dev raw[MV_MAX_DEVICES];
    int n = collect_devices(raw, MV_MAX_DEVICES);
    if (!buf)
        return n;
    if (n > cap)
        n = cap;
    for (int i = 0; i < n; i++) {
        mv_device_info *d = &buf[i];
        memset(d, 0, sizeof(*d));
        clGetPlatformInfo(raw[i].plat, CL_PLATFORM_NAME, sizeof d->platform, d->platform, NULL);
        clGetDeviceInfo(raw[i].dev, CL_DEVICE_NAME, sizeof d->name, d->name, NULL);
        cl_uint cu = 0, mhz = 0;
        cl_ulong gm = 0;
        cl_device_type ty = 0;
        clGetDeviceInfo(raw[i].dev, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof cu, &cu, NULL);
        clGetDeviceInfo(raw[i].dev, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof mhz, &mhz, NULL);
        clGetDeviceInfo(raw[i].dev, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof gm, &gm, NULL);
        clGetDeviceInfo(raw[i].dev, CL_DEVICE_TYPE, sizeof ty, &ty, NULL);
        d->compute_units = cu;
        d->clock_mhz = mhz;
        d->global_mem = (unsigned long long)gm;
        d->is_gpu = (ty & CL_DEVICE_TYPE_GPU) ? 1 : 0;
    }
    return n;
}

static void mv_grind_opts_default(mv_grind_opts *o)
{
    memset(o, 0, sizeof(*o));
    o->device_index = -1;
    o->iters_per_thread = 512;
    o->batch = 64;
}

static void *read_entire(const char *path, size_t *len)
{
    FILE *f = fopen(path, "rb");
    if (!f)
        return NULL;
    fseek(f, 0, SEEK_END);
    long n = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (n < 0) {
        fclose(f);
        return NULL;
    }
    char *b = (char *)malloc((size_t)n + 1);
    if (!b || fread(b, 1, (size_t)n, f) != (size_t)n) {
        free(b);
        fclose(f);
        return NULL;
    }
    fclose(f);
    b[n] = 0;
    if (len)
        *len = (size_t)n;
    return b;
}

static char *read_kernel_file(const char *dir, const char *name)
{
    char path[1024];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    return (char *)read_entire(path, NULL);
}

#define MV_NPARTS 4

static char *build_source(const mv_grind_opts *o)
{
    static const char *files[MV_NPARTS] = {
        "kernels/mv_prelude.h",
        "build/monocypher_x25519.inc",
        "kernels/mv_edwards.h",
        "kernels/mv_grind.cl",
    };
    const char *parts[MV_NPARTS];
    char *owned[MV_NPARTS] = {NULL, NULL, NULL, NULL};

    if (o->kernel_dir) {
        for (int i = 0; i < MV_NPARTS; i++) {
            owned[i] = read_kernel_file(o->kernel_dir, files[i]);
            if (!owned[i]) {
                for (int j = 0; j < MV_NPARTS; j++)
                    free(owned[j]);
                return NULL;
            }
            parts[i] = owned[i];
        }
    } else {
        parts[0] = MV_SRC_PRELUDE;
        parts[1] = MV_SRC_X25519;
        parts[2] = MV_SRC_EDWARDS;
        parts[3] = MV_SRC_GRIND;
    }

    size_t total = 0;
    for (int i = 0; i < MV_NPARTS; i++)
        total += strlen(parts[i]) + 2;
    char *src = (char *)malloc(total + 1);
    if (src) {
        src[0] = 0;
        for (int i = 0; i < MV_NPARTS; i++) {
            strcat(src, parts[i]);
            strcat(src, "\n");
        }
    }
    for (int i = 0; i < MV_NPARTS; i++)
        free(owned[i]);
    return src;
}

static uint64_t fnv1a(const char *s, uint64_t h)
{
    while (*s) {
        h ^= (unsigned char)*s++;
        h *= 1099511628211ULL;
    }
    return h;
}

static void cache_path(char *out, size_t outlen, cl_device_id id, const char *src, const char *opts)
{
    const char *dir = getenv("MV_CACHE_DIR");
    char fallback[512];
    if (!dir || !*dir) {
        const char *home = getenv("HOME");
        if (home && *home)
            snprintf(fallback, sizeof fallback, "%s/.cache/meshvanity", home);
        else
            snprintf(fallback, sizeof fallback, "/var/cache/meshvanity");
        dir = fallback;
    }

    char devname[256] = {0}, devver[256] = {0}, drvver[256] = {0};
    clGetDeviceInfo(id, CL_DEVICE_NAME, sizeof devname, devname, NULL);
    clGetDeviceInfo(id, CL_DEVICE_VERSION, sizeof devver, devver, NULL);
    clGetDeviceInfo(id, CL_DRIVER_VERSION, sizeof drvver, drvver, NULL);

    uint64_t h = 1469598103934665603ULL;
    h = fnv1a(src, h);
    h = fnv1a(opts, h);
    h = fnv1a(devname, h);
    h = fnv1a(devver, h);
    h = fnv1a(drvver, h);

    char tmp[512];
    snprintf(tmp, sizeof tmp, "%s", dir);
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;
            mkdir(tmp, 0700);
            *p = '/';
        }
    }
    mkdir(tmp, 0700);

    snprintf(out, outlen, "%s/prog-%016llx.bin", dir, (unsigned long long)h);
}

static cl_program program_from_cache(cl_context ctx, cl_device_id id, const char *path)
{
    size_t len = 0;
    unsigned char *buf = (unsigned char *)read_entire(path, &len);
    if (!buf || len == 0) {
        free(buf);
        return NULL;
    }

    const unsigned char *ptr = buf;
    cl_int bin_status = 0, e = 0;
    cl_program p = clCreateProgramWithBinary(ctx, 1, &id, &len, &ptr, &bin_status, &e);
    free(buf);
    if (e != CL_SUCCESS || bin_status != CL_SUCCESS) {
        if (p)
            clReleaseProgram(p);
        return NULL;
    }
    return p;
}

static void program_to_cache(cl_program prog, const char *path)
{
    size_t nbin = 0;
    if (clGetProgramInfo(prog, CL_PROGRAM_BINARY_SIZES, sizeof nbin, &nbin, NULL) != CL_SUCCESS || nbin == 0)
        return;

    unsigned char *buf = (unsigned char *)malloc(nbin);
    if (!buf)
        return;
    unsigned char *ptrs[1] = {buf};
    if (clGetProgramInfo(prog, CL_PROGRAM_BINARIES, sizeof ptrs, ptrs, NULL) != CL_SUCCESS) {
        free(buf);
        return;
    }

    char tmp[700];
    int n = snprintf(tmp, sizeof tmp, "%s.tmp%d", path, (int)getpid());
    if (n < 0 || (size_t)n >= sizeof tmp) {
        free(buf);
        return;
    }
    FILE *f = fopen(tmp, "wb");
    if (f) {
        size_t wrote = fwrite(buf, 1, nbin, f);
        fclose(f);
        if (wrote == nbin)
            rename(tmp, path);
        else
            remove(tmp);
    }
    free(buf);
}

static int dev_init(struct mv_dev *d, cl_device_id id, const char *src, const mv_grind_opts *o, char *err, size_t errlen)
{
    cl_int e;
    d->id = id;
    clGetDeviceInfo(id, CL_DEVICE_NAME, sizeof d->name, d->name, NULL);

    d->ctx = clCreateContext(NULL, 1, &id, NULL, NULL, &e);
    CHECK(e, "clCreateContext");
    d->q = clCreateCommandQueue(d->ctx, id, 0, &e);
    CHECK(e, "clCreateCommandQueue");

    char buildopts[128];
    snprintf(buildopts, sizeof buildopts, "-cl-std=CL1.2 -cl-mad-enable -D MV_BATCH=%u", o->batch ? o->batch : 1u);

    char cpath[600];
    int cached = 0;
    cache_path(cpath, sizeof cpath, id, src, buildopts);
    d->prog = program_from_cache(d->ctx, id, cpath);
    if (d->prog) {
        cached = 1;
    } else {
        d->prog = clCreateProgramWithSource(d->ctx, 1, &src, NULL, &e);
        CHECK(e, "clCreateProgramWithSource");
    }

    e = clBuildProgram(d->prog, 1, &id, buildopts, NULL, NULL);
    if (e != CL_SUCCESS && cached) {
        clReleaseProgram(d->prog);
        remove(cpath);
        cached = 0;
        d->prog = clCreateProgramWithSource(d->ctx, 1, &src, NULL, &e);
        CHECK(e, "clCreateProgramWithSource");
        e = clBuildProgram(d->prog, 1, &id, buildopts, NULL, NULL);
    }
    if (e != CL_SUCCESS) {
        size_t logn = 0;
        clGetProgramBuildInfo(d->prog, id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logn);
        char *log = (char *)malloc(logn + 1);
        if (log) {
            clGetProgramBuildInfo(d->prog, id, CL_PROGRAM_BUILD_LOG, logn, log, NULL);
            log[logn] = 0;
            snprintf(err, errlen, "build failed on %s:\n%s", d->name, log);
            free(log);
        } else {
            snprintf(err, errlen, "build failed on %s (no log)", d->name);
        }
        return -1;
    }

    if (!cached)
        program_to_cache(d->prog, cpath);

    d->k_grind = clCreateKernel(d->prog, "mv_grind", &e);
    CHECK(e, "clCreateKernel(mv_grind)");
    d->k_batch = clCreateKernel(d->prog, "mv_grind_batch", &e);
    CHECK(e, "clCreateKernel(mv_grind_batch)");
    d->k_derive = clCreateKernel(d->prog, "mv_derive", &e);
    CHECK(e, "clCreateKernel(mv_derive)");

    cl_uint cu = 1;
    clGetDeviceInfo(id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof cu, &cu, NULL);
    size_t wgmul = 64;
    clGetKernelWorkGroupInfo(d->k_grind, id, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof wgmul, &wgmul, NULL);
    if (wgmul == 0)
        wgmul = 64;
    d->global_size = o->global_size ? o->global_size : (size_t)cu * wgmul * 32;
    d->ramp = 1024;
    d->iters = o->iters_per_thread ? o->iters_per_thread : 16;

    unsigned b = o->batch ? o->batch : 1u;
    if (d->iters % b)
        d->iters += b - (d->iters % b);

    d->n_jumps = 1;
    while (((size_t)1 << d->n_jumps) < d->global_size && d->n_jumps < 40)
        d->n_jumps++;

    d->m_base = clCreateBuffer(d->ctx, CL_MEM_READ_ONLY, 32, NULL, &e);
    CHECK(e, "alloc base");
    d->m_targets = clCreateBuffer(d->ctx, CL_MEM_READ_ONLY, sizeof(cl_uint) * MV_MAX_TARGETS, NULL, &e);
    CHECK(e, "alloc targets");
    d->m_count = clCreateBuffer(d->ctx, CL_MEM_READ_WRITE, sizeof(cl_uint), NULL, &e);
    CHECK(e, "alloc count");
    d->m_k = clCreateBuffer(d->ctx, CL_MEM_WRITE_ONLY, sizeof(cl_ulong) * MV_OUT_CAP, NULL, &e);
    CHECK(e, "alloc k");
    d->m_id = clCreateBuffer(d->ctx, CL_MEM_WRITE_ONLY, sizeof(cl_uint) * MV_OUT_CAP, NULL, &e);
    CHECK(e, "alloc id");

    d->m_pstart = clCreateBuffer(d->ctx, CL_MEM_READ_ONLY, sizeof(cl_int) * MV_GE_INTS, NULL, &e);
    CHECK(e, "alloc pstart");
    d->m_jumps = clCreateBuffer(d->ctx, CL_MEM_READ_ONLY, sizeof(cl_int) * MV_GE_INTS * d->n_jumps, NULL, &e);
    CHECK(e, "alloc jumps");
    d->m_step = clCreateBuffer(d->ctx, CL_MEM_READ_ONLY, sizeof(cl_int) * MV_GE_INTS, NULL, &e);
    CHECK(e, "alloc step");
    return 0;
}

static void dev_free(struct mv_dev *d)
{
    if (d->m_base)
        clReleaseMemObject(d->m_base);
    if (d->m_targets)
        clReleaseMemObject(d->m_targets);
    if (d->m_count)
        clReleaseMemObject(d->m_count);
    if (d->m_k)
        clReleaseMemObject(d->m_k);
    if (d->m_id)
        clReleaseMemObject(d->m_id);
    if (d->m_pstart)
        clReleaseMemObject(d->m_pstart);
    if (d->m_jumps)
        clReleaseMemObject(d->m_jumps);
    if (d->m_step)
        clReleaseMemObject(d->m_step);
    if (d->k_grind)
        clReleaseKernel(d->k_grind);
    if (d->k_batch)
        clReleaseKernel(d->k_batch);
    if (d->k_derive)
        clReleaseKernel(d->k_derive);
    if (d->prog)
        clReleaseProgram(d->prog);
    if (d->q)
        clReleaseCommandQueue(d->q);
    if (d->ctx)
        clReleaseContext(d->ctx);
    memset(d, 0, sizeof(*d));
}

static void mv_close(mv_ctx *c)
{
    if (!c)
        return;
    for (int i = 0; i < c->n_dev; i++)
        dev_free(&c->dev[i]);
    free(c);
}

static mv_ctx *mv_open(const mv_grind_opts *opts, char *err, size_t errlen)
{
    struct raw_dev raw[MV_MAX_DEVICES];
    int nraw = collect_devices(raw, MV_MAX_DEVICES);
    if (nraw == 0) {
        snprintf(err, errlen, "no OpenCL devices found (is an ICD installed?)");
        return NULL;
    }

    mv_ctx *c = (mv_ctx *)calloc(1, sizeof(mv_ctx));
    if (!c) {
        snprintf(err, errlen, "out of memory");
        return NULL;
    }
    c->opts = *opts;

    char *src = build_source(opts);
    if (!src) {
        snprintf(err, errlen, "could not assemble kernel source%s", opts->kernel_dir ? " (check --kernel-dir)" : "");
        free(c);
        return NULL;
    }

    for (int i = 0; i < nraw; i++) {
        if (opts->device_index >= 0 && i != opts->device_index)
            continue;
        if (opts->device_index < 0) {
            cl_device_type ty = 0;
            clGetDeviceInfo(raw[i].dev, CL_DEVICE_TYPE, sizeof ty, &ty, NULL);
            if (!(ty & CL_DEVICE_TYPE_GPU))
                continue;
        }
        struct mv_dev *d = &c->dev[c->n_dev];
        if (dev_init(d, raw[i].dev, src, opts, err, errlen) != 0) {
            dev_free(d);
            free(src);
            mv_close(c);
            return NULL;
        }
        d->index = i;
        c->n_dev++;
    }
    free(src);

    if (c->n_dev == 0) {
        snprintf(err, errlen, "no usable device (index %d)", opts->device_index);
        free(c);
        return NULL;
    }
    return c;
}

static int mv_selftest_gpu(mv_ctx *c, unsigned n, char *err, size_t errlen)
{
    uint8_t *sk = (uint8_t *)malloc((size_t)n * 32);
    uint8_t *pk_gpu = (uint8_t *)malloc((size_t)n * 32);
    uint32_t *id_gpu = (uint32_t *)malloc((size_t)n * 4);
    if (!sk || !pk_gpu || !id_gpu) {
        free(sk);
        free(pk_gpu);
        free(id_gpu);
        snprintf(err, errlen, "out of memory");
        return -1;
    }

    uint64_t x = 0x9E3779B97F4A7C15ull;
    for (unsigned i = 0; i < n * 32; i++) {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        sk[i] = (uint8_t)(x >> 24);
    }
    for (unsigned i = 0; i < n; i++)
        mv_clamp(sk + (size_t)i * 32);

    for (int di = 0; di < c->n_dev; di++) {
        struct mv_dev *d = &c->dev[di];
        cl_int e;
        cl_mem m_sk = clCreateBuffer(d->ctx, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, (size_t)n * 32, sk, &e);
        CHECK(e, "selftest alloc sk");
        cl_mem m_pk = clCreateBuffer(d->ctx, CL_MEM_WRITE_ONLY, (size_t)n * 32, NULL, &e);
        CHECK(e, "selftest alloc pk");
        cl_mem m_id = clCreateBuffer(d->ctx, CL_MEM_WRITE_ONLY, (size_t)n * 4, NULL, &e);
        CHECK(e, "selftest alloc id");

        cl_uint cn = n;
        clSetKernelArg(d->k_derive, 0, sizeof(cl_mem), &m_sk);
        clSetKernelArg(d->k_derive, 1, sizeof(cl_mem), &m_pk);
        clSetKernelArg(d->k_derive, 2, sizeof(cl_mem), &m_id);
        clSetKernelArg(d->k_derive, 3, sizeof(cl_uint), &cn);

        size_t gs = n;
        e = clEnqueueNDRangeKernel(d->q, d->k_derive, 1, NULL, &gs, NULL, 0, NULL, NULL);
        CHECK(e, "selftest enqueue");
        e = clEnqueueReadBuffer(d->q, m_pk, CL_TRUE, 0, (size_t)n * 32, pk_gpu, 0, NULL, NULL);
        CHECK(e, "selftest read pk");
        e = clEnqueueReadBuffer(d->q, m_id, CL_TRUE, 0, (size_t)n * 4, id_gpu, 0, NULL, NULL);
        CHECK(e, "selftest read id");
        clFinish(d->q);
        clReleaseMemObject(m_sk);
        clReleaseMemObject(m_pk);
        clReleaseMemObject(m_id);

        for (unsigned i = 0; i < n; i++) {
            uint8_t want[32];
            crypto_x25519_public_key(want, sk + (size_t)i * 32);
            if (memcmp(want, pk_gpu + (size_t)i * 32, 32) != 0) {
                snprintf(err, errlen, "device %d (%s): pubkey mismatch at key %u", di, d->name, i);
                free(sk);
                free(pk_gpu);
                free(id_gpu);
                return -1;
            }
            if (mv_crc32_core(want, 32) != id_gpu[i]) {
                snprintf(err, errlen, "device %d (%s): nodenum mismatch at key %u", di, d->name, i);
                free(sk);
                free(pk_gpu);
                free(id_gpu);
                return -1;
            }
        }
    }
    free(sk);
    free(pk_gpu);
    free(id_gpu);
    return 0;
}

/* Key security rests on this read: a low-entropy seed (time, PID, counter)
 * makes every key this grinder produces re-derivable by anyone who brute-forces
 * the seed. Only a real CSPRNG here: getrandom() blocks until the pool is
 * seeded; both paths fail closed. */
static int fill_random(uint8_t *buf, size_t n, char *err, size_t errlen)
{
    int have = 0;

#ifdef MV_HAVE_GETRANDOM
    size_t off = 0;
    while (off < n) {
        ssize_t r = getrandom(buf + off, n - off, 0);
        if (r < 0) {
            if (errno == EINTR)
                continue;
            break;
        }
        off += (size_t)r;
    }
    have = (off == n);
#endif

    if (!have) {
        FILE *f = fopen("/dev/urandom", "rb");
        if (!f) {
            snprintf(err, errlen,
                     "no CSPRNG available (getrandom failed and "
                     "/dev/urandom will not open)");
            return -1;
        }
        size_t got = fread(buf, 1, n, f);
        fclose(f);
        if (got != n) {
            snprintf(err, errlen, "short read from /dev/urandom");
            return -1;
        }
    }

    int nonzero = 0;
    for (size_t i = 0; i < n; i++)
        nonzero |= buf[i];
    if (!nonzero) {
        snprintf(err, errlen, "CSPRNG returned all zeros, refusing to generate a key");
        return -1;
    }
    return 0;
}

static long mv_run(mv_ctx *c, const mv_pattern *pat, mv_hit_cb on_hit, mv_progress_cb on_prog, void *user, volatile int *stop,
                   char *err, size_t errlen)
{
    if (pat->n_targets > MV_MAX_TARGETS) {
        snprintf(err, errlen, "too many targets (%u > %d)", pat->n_targets, MV_MAX_TARGETS);
        return -1;
    }

    uint8_t base[MV_MAX_DEVICES][32];
    for (int i = 0; i < c->n_dev; i++) {
        if (fill_random(base[i], 32, err, errlen) != 0)
            return -1;
        mv_clamp(base[i]);
    }

    for (int i = 0; i < c->n_dev; i++) {
        struct mv_dev *d = &c->dev[i];
        cl_uint zero = 0;
        cl_int e;
        e = clEnqueueWriteBuffer(d->q, d->m_base, CL_TRUE, 0, 32, base[i], 0, NULL, NULL);
        CHECK(e, "write base");
        e = clEnqueueWriteBuffer(d->q, d->m_targets, CL_TRUE, 0, sizeof(cl_uint) * pat->n_targets, pat->targets, 0, NULL, NULL);
        CHECK(e, "write targets");
        e = clEnqueueWriteBuffer(d->q, d->m_count, CL_TRUE, 0, sizeof(cl_uint), &zero, 0, NULL, NULL);
        CHECK(e, "write count");

        cl_uint n_t = pat->n_targets, msk = pat->mask, cap = MV_OUT_CAP, it = d->iters;
        clSetKernelArg(d->k_grind, 0, sizeof(cl_mem), &d->m_base);
        clSetKernelArg(d->k_grind, 2, sizeof(cl_uint), &it);
        clSetKernelArg(d->k_grind, 3, sizeof(cl_mem), &d->m_targets);
        clSetKernelArg(d->k_grind, 4, sizeof(cl_uint), &n_t);
        clSetKernelArg(d->k_grind, 5, sizeof(cl_uint), &msk);
        clSetKernelArg(d->k_grind, 6, sizeof(cl_mem), &d->m_count);
        clSetKernelArg(d->k_grind, 7, sizeof(cl_mem), &d->m_k);
        clSetKernelArg(d->k_grind, 8, sizeof(cl_mem), &d->m_id);
        clSetKernelArg(d->k_grind, 9, sizeof(cl_uint), &cap);

        cl_uint nj = d->n_jumps;
        clSetKernelArg(d->k_batch, 0, sizeof(cl_mem), &d->m_pstart);
        clSetKernelArg(d->k_batch, 1, sizeof(cl_mem), &d->m_jumps);
        clSetKernelArg(d->k_batch, 2, sizeof(cl_uint), &nj);
        clSetKernelArg(d->k_batch, 3, sizeof(cl_mem), &d->m_step);
        clSetKernelArg(d->k_batch, 5, sizeof(cl_uint), &it);
        clSetKernelArg(d->k_batch, 6, sizeof(cl_mem), &d->m_targets);
        clSetKernelArg(d->k_batch, 7, sizeof(cl_uint), &n_t);
        clSetKernelArg(d->k_batch, 8, sizeof(cl_uint), &msk);
        clSetKernelArg(d->k_batch, 9, sizeof(cl_mem), &d->m_count);
        clSetKernelArg(d->k_batch, 10, sizeof(cl_mem), &d->m_k);
        clSetKernelArg(d->k_batch, 11, sizeof(cl_mem), &d->m_id);
        clSetKernelArg(d->k_batch, 12, sizeof(cl_uint), &cap);
    }

    uint64_t k_next[MV_MAX_DEVICES] = {0};
    uint64_t tried = 0;
    long hits = 0;
    const double t0 = now_sec();
    double t_report = t0;

    while (!(stop && *stop)) {
        for (int i = 0; i < c->n_dev; i++) {
            struct mv_dev *d = &c->dev[i];
            cl_ulong kb = (cl_ulong)k_next[i];
            cl_kernel kern;

            if (c->opts.use_ladder) {
                clSetKernelArg(d->k_grind, 1, sizeof(cl_ulong), &kb);
                kern = d->k_grind;
            } else {
                int32_t p_start[MV_GE_INTS], step[MV_GE_INTS];
                int32_t *jumps = (int32_t *)malloc(sizeof(int32_t) * MV_GE_INTS * d->n_jumps);
                if (!jumps) {
                    snprintf(err, errlen, "out of memory building walk");
                    return -1;
                }
                if (mv_build_walk(p_start, jumps, d->n_jumps, step, base[i], k_next[i], d->iters) != 0) {
                    free(jumps);
                    snprintf(err, errlen, "mv_build_walk failed");
                    return -1;
                }
                cl_int we;
                we = clEnqueueWriteBuffer(d->q, d->m_pstart, CL_TRUE, 0, sizeof(int32_t) * MV_GE_INTS, p_start, 0, NULL, NULL);
                if (we == CL_SUCCESS)
                    we = clEnqueueWriteBuffer(d->q, d->m_jumps, CL_TRUE, 0, sizeof(int32_t) * MV_GE_INTS * d->n_jumps, jumps, 0,
                                              NULL, NULL);
                if (we == CL_SUCCESS)
                    we = clEnqueueWriteBuffer(d->q, d->m_step, CL_TRUE, 0, sizeof(int32_t) * MV_GE_INTS, step, 0, NULL, NULL);
                free(jumps);
                CHECK(we, "upload walk");
                kern = d->k_batch;
                clSetKernelArg(kern, 4, sizeof(cl_ulong), &kb);
            }

            const size_t gs = d->ramp < d->global_size ? d->ramp : d->global_size;
            cl_int e = clEnqueueNDRangeKernel(d->q, kern, 1, NULL, &gs, NULL, 0, NULL, NULL);
            CHECK(e, "enqueue grind");
            if (d->ramp < d->global_size)
                d->ramp *= 4;
            k_next[i] += (uint64_t)gs * d->iters;
            tried += (uint64_t)gs * d->iters;
        }

        for (int i = 0; i < c->n_dev; i++) {
            struct mv_dev *d = &c->dev[i];
            cl_int e = clFinish(d->q);
            CHECK(e, "clFinish");

            cl_uint count = 0;
            e = clEnqueueReadBuffer(d->q, d->m_count, CL_TRUE, 0, sizeof count, &count, 0, NULL, NULL);
            CHECK(e, "read count");
            if (count == 0)
                continue;

            cl_uint n_read = count > MV_OUT_CAP ? MV_OUT_CAP : count;
            cl_ulong ks[MV_OUT_CAP];
            cl_uint ids[MV_OUT_CAP];
            e = clEnqueueReadBuffer(d->q, d->m_k, CL_TRUE, 0, sizeof(cl_ulong) * n_read, ks, 0, NULL, NULL);
            CHECK(e, "read k");
            e = clEnqueueReadBuffer(d->q, d->m_id, CL_TRUE, 0, sizeof(cl_uint) * n_read, ids, 0, NULL, NULL);
            CHECK(e, "read ids");
            cl_uint zero = 0;
            clEnqueueWriteBuffer(d->q, d->m_count, CL_TRUE, 0, sizeof zero, &zero, 0, NULL, NULL);

            for (cl_uint h = 0; h < n_read; h++) {
                mv_hit hit;
                memset(&hit, 0, sizeof hit);
                mv_secret_at(hit.secret, base[i], (uint64_t)ks[h]);
                crypto_x25519_public_key(hit.pub, hit.secret);
                hit.nodenum = mv_crc32_core(hit.pub, 32);
                hit.device_index = d->index;

                if (hit.nodenum != ids[h] || !mv_pattern_match(pat, hit.nodenum)) {
                    snprintf(err, errlen,
                             "device %d (%s) reported id %08x but CPU re-derivation gives %08x: "
                             "kernel is wrong, refusing the result",
                             d->index, d->name, (unsigned)ids[h], hit.nodenum);
                    return -1;
                }
                if (!mv_is_clamped(hit.secret)) {
                    snprintf(err, errlen, "produced an unclamped secret, refusing");
                    return -1;
                }

                hits++;
                if (on_hit && on_hit(&hit, user))
                    return hits;
            }
        }

        const double t = now_sec();
        if (on_prog && t - t_report >= 1.0) {
            on_prog(tried, (double)tried / (t - t0), user);
            t_report = t;
        }
    }
    return hits;
}

/* ------------------------------------------------------------------------- */
/* CLI                                                                        */
/* ------------------------------------------------------------------------- */

static volatile int g_stop = 0;
static void on_sigint(int sig)
{
    (void)sig;
    g_stop = 1;
}

struct app {
    const char *outpath;
    int found_target;
    long want;
};

static void b64(char *dst, const uint8_t *src, size_t n)
{
    static const char *T = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    size_t o = 0;
    for (size_t i = 0; i < n; i += 3) {
        unsigned v = src[i] << 16;
        if (i + 1 < n)
            v |= src[i + 1] << 8;
        if (i + 2 < n)
            v |= src[i + 2];
        dst[o++] = T[(v >> 18) & 63];
        dst[o++] = T[(v >> 12) & 63];
        dst[o++] = (i + 1 < n) ? T[(v >> 6) & 63] : '=';
        dst[o++] = (i + 2 < n) ? T[v & 63] : '=';
    }
    dst[o] = 0;
}

static void hexify(char *dst, const uint8_t *src, size_t n)
{
    static const char *H = "0123456789abcdef";
    for (size_t i = 0; i < n; i++) {
        dst[2 * i] = H[src[i] >> 4];
        dst[2 * i + 1] = H[src[i] & 15];
    }
    dst[2 * n] = 0;
}

static int on_hit(const mv_hit *h, void *user)
{
    struct app *a = (struct app *)user;
    char skh[65], pkh[65], skb[64], pkb[64], idb[10];
    hexify(skh, h->secret, 32);
    hexify(pkh, h->pub, 32);
    b64(skb, h->secret, 32);
    b64(pkb, h->pub, 32);
    snprintf(idb, sizeof idb, "!%08x", h->nodenum);

    const mv_color hc = mv_color_of_nodenum(h->nodenum);
    char cpretty[64], cplain[8];
    mv_color_fmt(cpretty, sizeof cpretty, hc, isatty(fileno(stdout)));
    mv_color_fmt(cplain, sizeof cplain, hc, 0);

    printf("\n\n=== HIT (device %d) ===\n", h->device_index);
    printf("node id     : %s\n", idb);
    printf("app color   : %s\n", cpretty);
    printf("private key : %s\n", skh);
    printf("public key  : %s\n", pkh);
    printf("priv base64 : %s\n", skb);
    printf("pub  base64 : %s\n", pkb);
    printf("clamped     : %s\n", mv_is_clamped(h->secret) ? "yes" : "NO");
    fflush(stdout);

    if (a->outpath) {
        FILE *f = fopen(a->outpath, "a");
        if (f) {
            fprintf(f,
                    "node_id=%s\napp_color=%s\nprivate_key_hex=%s\npublic_key_hex=%s\n"
                    "private_key_b64=%s\npublic_key_b64=%s\n\n",
                    idb, cplain, skh, pkh, skb, pkb);
            fclose(f);
        }
    }
    a->found_target++;
    return (a->want > 0 && a->found_target >= a->want) ? 1 : 0;
}

static double g_expected = 0;
static double g_bench_secs = 0;
static double g_bench_rate = 0;

static void on_prog(uint64_t tried, double rate, void *user)
{
    (void)user;
    double eta = rate > 0 ? (g_expected - (double)tried) / rate : 0;
    if (eta < 0)
        eta = 0;
    printf("\r%" PRIu64 " keys  %.2f M/s  eta ~%.1f min      ", tried, rate / 1e6, eta / 60.0);
    fflush(stdout);

    if (g_bench_secs > 0 && rate > 0) {
        g_bench_rate = rate;
        if ((double)tried / rate >= g_bench_secs)
            g_stop = 1;
    }
}

static void usage(const char *p)
{
    fprintf(stderr,
            "mvgrind: GPU vanity NodeNum grinder for Meshtastic\n\n"
            "usage: %s [pattern] [options]\n\n"
            "pattern:\n"
            "  dc801051         exact node id (a leading '!' is fine)\n"
            "  dc80             prefix: top 4 hex digits\n"
            "  dc80****         explicit wildcard nibbles ('*', '?', '.')\n"
            "  dc80,801f,d0f0   a set; any match wins, at no extra cost\n"
            "  (omit it entirely when you only care about --color)\n\n"
            "color:\n"
            "  the apps paint a node with the low 24 bits of its id as RGB, so\n"
            "  --color is a pattern over those bits and combines with the one above\n"
            "  (id nibbles 3-8 ARE the channels: 'dc80' already pins red to 0x80)\n\n"
            "options:\n"
            "  -c, --color SPEC   #rrggbb, #rgb, or any CSS color name\n"
            "      --tol N        accept +/-N per channel (default 0, exact)\n"
            "  -d, --device N     device index (default: every GPU)\n"
            "  -n, --count N      stop after N hits (default 1, 0 = never)\n"
            "  -o, --out FILE     append hits to FILE (default found.txt)\n"
            "  -i, --iters N      candidates per work-item per dispatch (default 16)\n"
            "  -g, --global N     work-items per dispatch (default: auto)\n"
            "  -k, --kernel-dir D load kernels from D instead of the built-in copies\n"
            "  -l, --list         list OpenCL devices and exit\n"
            "  -t, --selftest     run CPU vectors + CPU/GPU differential and exit\n"
            "  -b, --bench SECS   measure throughput for SECS and exit\n"
            "      --batch N      batched-inversion width (default 8, 1 = off)\n"
            "      --ladder       use the old per-candidate ladder kernel\n"
            "  -h, --help         this text\n",
            p);
}

static const char *need_value(const char *opt, const char *val)
{
    if (!val) {
        fprintf(stderr, "%s needs a value\n", opt);
        exit(2);
    }
    return val;
}

int main(int argc, char **argv)
{
    char err[4096] = {0};
    mv_grind_opts o;
    mv_grind_opts_default(&o);

    struct app app;
    memset(&app, 0, sizeof app);
    app.outpath = "found.txt";
    app.want = 1;

    const char *pattern = NULL;
    const char *colorspec = NULL;
    int tol = 0, have_tol = 0;
    int do_list = 0, do_selftest = 0;
    double bench = 0;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
#define NEXT() need_value(a, i + 1 < argc ? argv[++i] : NULL)
        if (!strcmp(a, "-l") || !strcmp(a, "--list"))
            do_list = 1;
        else if (!strcmp(a, "-t") || !strcmp(a, "--selftest"))
            do_selftest = 1;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage(argv[0]);
            return 0;
        } else if (!strcmp(a, "-d") || !strcmp(a, "--device"))
            o.device_index = atoi(NEXT());
        else if (!strcmp(a, "-i") || !strcmp(a, "--iters"))
            o.iters_per_thread = (unsigned)atoi(NEXT());
        else if (!strcmp(a, "-g") || !strcmp(a, "--global"))
            o.global_size = (size_t)atoll(NEXT());
        else if (!strcmp(a, "-k") || !strcmp(a, "--kernel-dir"))
            o.kernel_dir = NEXT();
        else if (!strcmp(a, "-n") || !strcmp(a, "--count"))
            app.want = atol(NEXT());
        else if (!strcmp(a, "-o") || !strcmp(a, "--out"))
            app.outpath = NEXT();
        else if (!strcmp(a, "-b") || !strcmp(a, "--bench"))
            bench = atof(NEXT());
        else if (!strcmp(a, "-c") || !strcmp(a, "--color"))
            colorspec = NEXT();
        else if (!strcmp(a, "--tol") || !strcmp(a, "--tolerance")) {
            tol = atoi(NEXT());
            have_tol = 1;
        } else if (!strcmp(a, "--ladder"))
            o.use_ladder = 1;
        else if (!strcmp(a, "--batch"))
            o.batch = (unsigned)atoi(NEXT());
        else if (a[0] == '-') {
            fprintf(stderr, "unknown option %s\n", a);
            return 1;
        } else
            pattern = a;
#undef NEXT
    }

    if (do_list) {
        mv_device_info info[16];
        int n = mv_list_devices(info, 16);
        if (n <= 0) {
            printf("no OpenCL devices found\n");
            return 1;
        }
        for (int i = 0; i < n; i++)
            printf("[%d] %s | %s | %u CU | %u MHz | %llu MB | %s\n", i, info[i].name, info[i].platform, info[i].compute_units,
                   info[i].clock_mhz, info[i].global_mem >> 20, info[i].is_gpu ? "GPU" : "CPU");
        return 0;
    }

    if (mv_selftest_cpu(err, sizeof err) != 0) {
        fprintf(stderr, "CPU selftest failed: %s\n", err);
        return 1;
    }

    if (do_selftest) {
        printf("cpu  : ok (RFC 7748 vectors, CRC-32 check value, clamp round-trip, color patterns)\n");
        mv_ctx *c = mv_open(&o, err, sizeof err);
        if (!c) {
            fprintf(stderr, "mv_open: %s\n", err);
            return 1;
        }
        if (mv_selftest_gpu(c, 8192, err, sizeof err) != 0) {
            fprintf(stderr, "GPU selftest FAILED: %s\n", err);
            mv_close(c);
            return 1;
        }
        printf("gpu  : ok (8192 random clamped keys, GPU == CPU pubkey and nodenum)\n");
        mv_close(c);
        return 0;
    }

    if (!pattern && !colorspec) {
        usage(argv[0]);
        return 1;
    }
    if (have_tol && !colorspec) {
        fprintf(stderr, "--tol only means something alongside --color\n");
        return 1;
    }

    mv_pattern pat;
    memset(&pat, 0, sizeof pat);
    char msg[256];
    if (pattern && mv_pattern_parse(&pat, pattern, msg, sizeof msg) != 0) {
        fprintf(stderr, "bad pattern: %s\n", msg);
        return 1;
    }

    mv_color color = {0, 0, 0};
    if (colorspec) {
        if (mv_color_parse(colorspec, &color, msg, sizeof msg) != 0) {
            fprintf(stderr, "bad color: %s\n", msg);
            mv_pattern_free(&pat);
            return 1;
        }
        if (mv_pattern_apply_color(&pat, color, tol, msg, sizeof msg) != 0) {
            fprintf(stderr, "%s\n", msg);
            mv_pattern_free(&pat);
            return 1;
        }
    }
    g_expected = mv_pattern_expected_keys(&pat);

    mv_ctx *c = mv_open(&o, err, sizeof err);
    if (!c) {
        fprintf(stderr, "mv_open: %s\n", err);
        mv_pattern_free(&pat);
        return 1;
    }

    if (mv_selftest_gpu(c, 2048, err, sizeof err) != 0) {
        fprintf(stderr, "GPU selftest FAILED, refusing to grind: %s\n", err);
        mv_close(c);
        mv_pattern_free(&pat);
        return 1;
    }

    mv_device_info info[16];
    int nd = mv_list_devices(info, 16);
    printf("pattern     : %s  (mask %08x, %u target%s)\n", pattern ? pattern : "(any id)", pat.mask, pat.n_targets,
           pat.n_targets == 1 ? "" : "s");
    if (colorspec) {
        char cb[64];
        mv_color_fmt(cb, sizeof cb, color, isatty(fileno(stdout)));
        printf("color       : %s", cb);
        if (tol)
            printf("   +/-%d per channel", tol);
        printf("\n");
    }
    printf("expected    : %.4g keys (mean)\n", g_expected);
    for (int i = 0; i < nd; i++)
        if (o.device_index < 0 ? info[i].is_gpu : i == o.device_index)
            printf("device [%d]  : %s (%u CU)\n", i, info[i].name, info[i].compute_units);
    printf("selftest    : passed\n\n");

    signal(SIGINT, on_sigint);
    if (bench > 0) {
        app.want = 0;
        g_bench_secs = bench;
        printf("benchmark mode: %.0fs\n", bench);
    }

    long hits = mv_run(c, &pat, on_hit, on_prog, &app, &g_stop, err, sizeof err);
    printf("\n");
    if (g_bench_secs > 0 && g_bench_rate > 0) {
        printf("benchmark: %.0f keys/sec (%.2f M/s)\n", g_bench_rate, g_bench_rate / 1e6);
        printf("expected time for a full 8-digit id: %.1f min (mean)\n", 4294967296.0 / g_bench_rate / 60.0);
    }
    if (hits < 0) {
        fprintf(stderr, "mv_run: %s\n", err);
        mv_close(c);
        mv_pattern_free(&pat);
        return 1;
    }
    printf("%ld hit%s\n", hits, hits == 1 ? "" : "s");
    mv_close(c);
    mv_pattern_free(&pat);
    return 0;
}
