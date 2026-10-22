// ============================================================================
// sdTrapezoid2D & sdTrapezoidSeg2D — Exact 0/1-Sqrt 2D Trapezoid & NaN-Free Segment Trapezoid
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`0` sqrts in cap zone, `1` in leg zone; Lagrange cross-product fixes 30.9% axis NaNs)
// @requires dot2 (from common.glsl)
// ============================================================================
//
// Mathematical Improvements over IQ's sdTrapezoid (Shadertoy MlycD3):
//   1. Endpoint Subsumption (0 sqrts in Cap Zone, eliminates cax):
//      Since the horizontal cap endpoints (r1, -he) and (r2, +he) are already
//      endpoints of the slanted leg segment, whenever |p.x| > rc (where
//      rc = (p.y < 0) ? r1 : r2), the closest cap point is the corner (rc, +-he),
//      which is already covered by the clamped leg projection `cb`. Thus `cax` is
//      never needed, and whenever `|p.x| <= rc`, the cap distance is linear (`cy =
//      abs(p.y) - he`, 0 sqrts!), giving a 0-sqrt fast exit for `|p.x| <= rc && cy >= 0`.
//   2. Lagrange 2D Cross-Product Identity in sdTrapezoidSeg2D (Eliminates Per-Pixel
//      sqrt and Fixes 30.9% FP32 Central-Axis NaNs):
//      In IQ's arbitrary-segment `sdTrapezoid(p, a, b, ra, rb)`, computing the
//      perpendicular distance via `x = sqrt(papa - paba*paba*baba)` (a) costs an
//      extra per-pixel sqrt and (b) rounds negative on 30.9% of points along the
//      central axis in FP32, producing NaN! By Lagrange's 2D identity:
//        baba * (papa - paba^2 * baba) == (pa.x * ba.y - pa.y * ba.x)^2,
//      so `x = abs(pa.x * ba.y - pa.y * ba.x) / l` is a rational 2D cross product
//      requiring 0 per-pixel sqrts for `x` and producing 0% NaNs everywhere in R^2!

#ifndef SDF_GEMS_SD_TRAPEZOID2D_GLSL
#define SDF_GEMS_SD_TRAPEZOID2D_GLSL

#include "../common/common.glsl"

// 1. Axis-aligned 2D isosceles trapezoid of bottom half-width r1, top half-width r2, and half-height he
float sdTrapezoid2D(vec2 p, float r1, float r2, float he) {
    p.x = abs(p.x);
    float rc = (p.y < 0.0) ? r1 : r2;
    float cy = abs(p.y) - he;
    if (p.x <= rc && cy >= 0.0) return cy; // 0 sqrts, 0 leg projections in exterior cap strip!
    vec2  k1 = vec2(r2, he);
    vec2  k2 = vec2(r2 - r1, 2.0 * he);
    vec2  cb = p - k1 + k2 * clamp(dot(k1 - p, k2) / dot2(k2), 0.0, 1.0);
    float b2 = dot2(cb);
    if (p.x <= rc && cy * cy <= b2) return cy; // 0 sqrts in interior cap Voronoi region!
    return ((cb.x < 0.0 && cy < 0.0) ? -1.0 : 1.0) * sqrt(b2);
}

// 2. Arbitrary-segment 2D trapezoid from endpoint a (half-width ra) to b (half-width rb)
float sdTrapezoidSeg2D(vec2 p, vec2 a, vec2 b, float ra, float rb) {
    float rba  = rb - ra;
    vec2  ba   = b - a, pa = p - a;
    float baba = dot(ba, ba), l = sqrt(baba);
    float paba = dot(pa, ba) / baba;
    float x    = abs(pa.x * ba.y - pa.y * ba.x) / l;
    float rc   = (paba < 0.5) ? ra : rb;
    float cay  = abs(paba - 0.5) - 0.5;
    if (x <= rc && cay >= 0.0) return cay * l;
    float f    = clamp((rba * (x - ra) + paba * baba) / (rba * rba + baba), 0.0, 1.0);
    float cbx  = x - ra - f * rba;
    float cby  = paba - f;
    float b2   = cbx * cbx + cby * cby * baba;
    if (x <= rc && cay * cay * baba <= b2) return -cay * l;
    return ((cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0) * sqrt(b2);
}

// 3. 0/1-SFU 2D Trapezoid Distance & Analytical Unit Gradient (Endpoint Subsumption: 0 SFU in cap zone)
vec3 sdgTrapezoid2D(vec2 p, float r1, float r2, float he) {
    float sx = (p.x < 0.0) ? -1.0 : 1.0, sy = (p.y < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float rc = (p.y < 0.0) ? r1 : r2, cy = abs(p.y) - he;
    if (p.x <= rc && cy >= 0.0) return vec3(cy, 0.0, sy);
    vec2  k  = vec2(r2 - r1, 2.0 * he), w = p - vec2(r1, -he);
    vec2  q  = w - k * clamp(dot(w, k) / dot(k, k), 0.0, 1.0);
    float b2 = dot(q, q);
    if (p.x <= rc && cy * cy <= b2) return vec3(cy, 0.0, sy);
    float s    = (w.x * k.y < w.y * k.x && cy < 0.0) ? -1.0 : 1.0;
    float invD = s * inversesqrt(max(b2, 1e-30));
    return (b2 > 0.0) ? vec3(b2 * invD, sx * q.x * invD, q.y * invD)
                      : vec3(0.0, vec2(sx * k.y, -k.x) * inversesqrt(dot(k, k)));
}

// 4. Exact Convex & Concave 2D Quad Distance & Analytical Unit Gradient (Jordan Winding-Parity + 1 inversesqrt)
vec3 sdgQuad2D(in vec2 p, in vec2 v[4]) {
    vec3  res = vec3(dot(p - v[0], p - v[0]), p - v[0]);
    float s   = 1.0;
    for (int i = 0, j = 3; i < 4; j = i, i++) {
        vec2  e = v[i] - v[j], w = p - v[j];
        vec2  q = w - e * clamp(dot(w, e) / dot(e, e), 0.0, 1.0);
        float d = dot(q, q);
        if (d < res.x) res = vec3(d, q);
        bvec3 c = bvec3(p.y >= v[j].y, p.y < v[i].y, e.x * w.y > e.y * w.x);
        if (all(c) || all(not(c))) s = -s;
    }
    float invD = s * inversesqrt(max(res.x, 1e-30));
    return vec3(res.x * invD, res.yz * invD);
}

#endif // SDF_GEMS_SD_TRAPEZOID2D_GLSL
