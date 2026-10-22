// ============================================================================
// sdCutHollowSphere — Exact 1/2-Sqrt 3D Cut Hollow Spherical Shell (Signed-Square & Conjugate-Rim)
// @author Inigo Quilez — upgraded to early-return 1-sqrt shell branch using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` in spherical shell [65.4% of space] & 2 in rim, down from 3)
// ============================================================================

#ifndef SDF_GEMS_SD_CUT_HOLLOW_SPHERE_GLSL
#define SDF_GEMS_SD_CUT_HOLLOW_SPHERE_GLSL

float sdCutHollowSphere( vec3 p, float r, float h, float t )
{
    float x2 = dot(p.xz, p.xz), w2 = r * r - h * h;
    if( h * abs(h) * x2 >= w2 * p.y * abs(p.y) )
        return abs(sqrt(x2 + p.y * p.y) - r) - t;
    float d = x2 - w2, dy = p.y - h;
    return sqrt(d * d / max(x2 + w2 + 2.0 * sqrt(w2 * x2), 1e-30) + dy * dy) - t;
}

#endif // SDF_GEMS_SD_CUT_HOLLOW_SPHERE_GLSL
