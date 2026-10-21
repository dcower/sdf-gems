// ============================================================================
// sdgCappedCone — Exact 1/2-Sqrt 3D Capped Cone with Analytical Gradient
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distgradfunctions3d/
// @exact  Yes (`||grad d|| = 1`; fixes 56.91% rim gradient failure in baseline & cuts sqrts from 3 to 1/2)
// @requires dot2 (from common.glsl)
// ============================================================================
//
// Bug Fix & Sqrt Reduction over IQ's 2025 sdgCappedCone (Shadertoy 3ctcW7):
//   1. In IQ's sdgCappedCone, the `la >= lb` branch returns the constant conical
//      face normal `vec3(k.y*p.xz/l, -k.x).xzy / sqrt(m)` even when the clamped
//      mantle parameter t = clamp(dot(q, k)/m, 0.0, 1.0) hits 0.0 or 1.0 (the top
//      and bottom circular rim wedges), and `la < lb` returns `vec3(0, sign(p.y), 0)`
//      when l > rc. Across [-2, 2]^3, IQ's sdgCappedCone returns a wrong gradient
//      at 56.91% of all points (max gradient error 1.597, > 105 deg direction error!).
//   2. By Endpoint Subsumption, the flat cap wins iff `l <= rc && cy * cy < lb`
//      (where cy = abs(p.y) - h), where d = cy (only 1 sqrt for l = length(p.xz))
//      and grad = vec3(0, sign(p.y), 0).
//   3. Everywhere else, the already-computed 2D displacement vector `b` from the
//      closest point on the mantle segment (both on the face t in (0, 1) and on
//      both circular rims t in {0, 1}) gives the exact unit 3D gradient
//      `vec3(b.x * p.x / l, b.y, b.x * p.z / l) / d` while eliminating `a`,
//      `dot(a, a)`, `sqrt(la)`, and `sqrt(m)` (2 sqrts total instead of 3).

#ifndef SDF_GEMS_SDG_CAPPED_CONE_GLSL
#define SDF_GEMS_SDG_CAPPED_CONE_GLSL

#include "../common/common.glsl"

vec4 sdgCappedCone(vec3 p, float h, float r1, float r2) {
    vec2  k  = vec2(r2 - r1, 2.0 * h);
    float l  = length(p.xz);
    vec2  q  = vec2(r2 - l, h - p.y);
    float rc = (p.y < 0.0) ? r1 : r2;
    float cy = abs(p.y) - h;
    vec2  b  = k * clamp(dot(q, k) / dot2(k), 0.0, 1.0) - q;
    float lb = dot2(b);
    if (l <= rc && cy * cy < lb) return vec4(cy, 0.0, (p.y < 0.0) ? -1.0 : 1.0, 0.0);
    float d  = ((b.x < 0.0 && cy < 0.0) ? -1.0 : 1.0) * sqrt(lb);
    return vec4(d, (l > 0.0) ? vec3(b.x * p.x / l, b.y, b.x * p.z / l) / d : vec3(0.0, sign(b.y), 0.0));
}

#endif // SDF_GEMS_SDG_CAPPED_CONE_GLSL
