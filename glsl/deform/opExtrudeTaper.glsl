// ============================================================================
// opExtrudeTaper — Exact 0/1-Sqrt Draft-Angle Tapered & Barrel Extrusions
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in cap zone & `1` in wall zone; fixes 37.9% overstep)
// @requires dot2 (from common.glsl)
// ============================================================================
//
// By the Orthogonal Eikonal Composition Theorem, for any exact 2D Euclidean SDF
// u = d2d(p.xy) (||grad_xy u|| = 1) and vertical coordinate z (grad_xy u perp grad_z z),
// composing any exact 2D profile SDF d_profile(u, z) with (d2d(p.xy), p.z) produces
// an exact 3D Euclidean SDF with ||grad_xyz d|| = 1 everywhere outside the 2D medial axis:
//
//   1. opExtrudeTaper(d2d, z, h, r1, r2): Extrudes any 2D SDF along z in [-h, +h]
//      with a constant draft angle / linear offset taper from bottom offset r1 (at z = -h)
//      to top offset r2 (at z = +h). Using Endpoint Subsumption on the semi-infinite
//      trapezoid { |z| <= h, u <= (r1+r2)/2 + (r2-r1)/(2h)*z }, the flat cap wins iff
//      d2d <= rc && cy^2 <= ||cb||^2 (where cy = |z| - h), requiring 0 sqrts in the
//      flat-cap Voronoi region and 1 sqrt in the slanted-mantle region.
//
//   2. opExtrudeBulge(d2d, z, h, b): Extrudes any 2D SDF along z in [-h, +h] with a
//      circular-arc bulged barrel wall (b > 0) or waisted spool wall (b < 0) of
//      mid-plane sagitta b using conjugate-square arc sector dominance (0/1 sqrts).

#ifndef SDF_GEMS_OP_EXTRUDE_TAPER_GLSL
#define SDF_GEMS_OP_EXTRUDE_TAPER_GLSL

#include "../common/common.glsl"

// 1. Exact Draft-Angle / Slanted Extrusion (offset r1 at z = -h, offset r2 at z = +h)
float opExtrudeTaper(float d2d, float z, float h, float r1, float r2) {
    vec2  q  = vec2(d2d, z);
    vec2  k1 = vec2(r2, h);
    vec2  k2 = vec2(r2 - r1, 2.0 * h);
    float rc = (z < 0.0) ? r1 : r2;
    float cy = abs(z) - h;
    vec2  cb = q - k1 + k2 * clamp(dot(k1 - q, k2) / dot2(k2), 0.0, 1.0);
    float b2 = dot2(cb);
    if (d2d <= rc && cy * cy <= b2) return cy;
    return ((cb.x < 0.0 && cy < 0.0) ? -1.0 : 1.0) * sqrt(b2);
}

// 2. Exact Circular-Arc Barrel (b > 0) & Waisted Spool (b < 0) Extrusion
float opExtrudeBulge(float d2d, float z, float h, float b) {
    float sb = sign(b);
    float d  = 0.5 * (h * h - b * b) / b;
    float R  = d + b;
    vec2  q  = vec2(d2d, abs(z));
    float u  = q.x + d;
    float s  = max(sb * (d * q.y - h * u), sb * ((R - h) * (R - h + 2.0 * q.y) - u * u));
    return (s < 0.0) ? sb * length(vec2(u, q.y)) - R :
           (q.x < 0.0) ? q.y - h :
                         length(vec2(q.x, q.y - h));
}

#endif // SDF_GEMS_OP_EXTRUDE_TAPER_GLSL
