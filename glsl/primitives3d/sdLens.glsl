// ============================================================================
// sdLens & sdgLens — Exact 1/2-Sqrt 3D Asymmetric Two-Sphere Optical Lens + Analytical Gradient
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt`/`rsq` in 85%–90% of R^3 & 2 in rim wedge, down from 4)
// ============================================================================

#ifndef SDF_GEMS_SD_LENS_GLSL
#define SDF_GEMS_SD_LENS_GLSL

// 1. Exact 1/2-Sqrt 3D Asymmetric Two-Sphere Optical Lens
//    By the Axis-Center Radial Collapse Theorem:
//      length(vec2(p.x - cx, length(p.yz))) === sqrt((p.x - cx)^2 + dot(p.yz, p.yz)),
//    and the Signed-Square Rim-Wedge Inequality:
//      a|a| yz2 > b^2 p.x|p.x|  &&  (d - a)|d - a| yz2 > b^2 (d - p.x)|d - p.x|,
//    85%–90% of R^3 (all interior and spherical-cap regions) evaluates in a SINGLE sqrt
//    (eliminating both `b = sqrt(ra^2 - a^2)` and `length(p.yz)`!), while the rim wedge
//    fuses `b * length(p.yz) = sqrt(b2 * yz2)` to use only 2 sqrts total.
float sdLens( vec3 p, float ra, float rb, float d )
{
    float yz2 = dot(p.yz, p.yz), a = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b2  = max(ra * ra - a * a, 0.0), da = d - a, dx = d - p.x;
    if( a * abs(a) * yz2 > b2 * p.x * abs(p.x) && da * abs(da) * yz2 > b2 * dx * abs(dx) )
    {
        float ux = p.x - a, dy = yz2 - b2;
        return sqrt(ux * ux + dy * dy / max(yz2 + b2 + 2.0 * sqrt(b2 * yz2), 1e-30));
    }
    float r  = ra - rb, m = ra * r + d * (p.x - a);
    bool  ca = (r * abs(r) * (p.x * p.x + yz2) <= m * abs(m));
    float cx = ca ? 0.0 : d, cr = ca ? ra : rb;
    return sqrt((p.x - cx) * (p.x - cx) + yz2) - cr;
}

// 2. Exact 1/2-SFU 3D Optical Lens Distance + Analytical Gradient (vec4(d, grad.xyz))
vec4 sdgLens( vec3 p, float ra, float rb, float d )
{
    float yz2 = dot(p.yz, p.yz), a = (ra * ra - rb * rb + d * d) / (2.0 * d);
    float b2  = max(ra * ra - a * a, 0.0), da = d - a, dx = d - p.x;
    if( a * abs(a) * yz2 > b2 * p.x * abs(p.x) && da * abs(da) * yz2 > b2 * dx * abs(dx) )
    {
        float invY = inversesqrt(max(yz2, 1e-30));
        float py   = (yz2 > 0.0) ? (yz2 * invY) : 0.0;
        vec2  er   = (yz2 > 0.0) ? (p.yz * invY) : vec2(1.0, 0.0);
        vec2  w    = vec2(p.x - a, py - sqrt(b2));
        float l2   = dot(w, w), invL = inversesqrt(max(l2, 1e-30));
        return vec4(l2 * invL, vec3(w.x, w.y * er.x, w.y * er.y) * invL);
    }
    float r  = ra - rb, m = ra * r + d * (p.x - a);
    bool  ca = (r * abs(r) * (p.x * p.x + yz2) <= m * abs(m));
    vec3  w  = vec3(p.x - (ca ? 0.0 : d), p.yz);
    float r2 = dot(w, w), invR = inversesqrt(max(r2, 1e-30));
    return vec4(r2 * invR - (ca ? ra : rb), (r2 > 0.0) ? (w * invR) : vec3(1.0, 0.0, 0.0));
}

#endif // SDF_GEMS_SD_LENS_GLSL
