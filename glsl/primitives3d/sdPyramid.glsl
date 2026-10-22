// ============================================================================
// sdPyramid — Exact 0/1-Sqrt 3D Square Pyramid (Three-Zone Euclidean Decomposition)
// @author Inigo Quilez — upgraded to 0/1-sqrt & base-bug-fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1` everywhere in R^3; 0 per-pixel sqrts in interior, base, & slanted faces)
// ============================================================================

#ifndef SDF_GEMS_SD_PYRAMID_GLSL
#define SDF_GEMS_SD_PYRAMID_GLSL

float sdPyramid( vec3 p, float h )
{
    vec2 u = abs(p.xz);
    u = (u.y > u.x) ? u.yx : u;
    u -= 0.5;

    float m2 = h * h + 0.25;
    float qy = h * p.y - 0.5 * u.x;
    float qz = h * u.x + 0.5 * p.y;

    if( (p.y <= 0.0 && u.x <= 0.0) || (qz <= 0.0 && p.y > 0.0 && m2 * p.y * p.y <= qz * qz) )
        return -p.y; // 0 sqrts in exterior & interior base-face prism!
    if( (qz <= 0.0 && p.y > 0.0) || (qy > 0.0 && -u.y * m2 - 0.5 * qy > 0.0) )
        return qz * inversesqrt(m2); // 0 per-pixel sqrts in slanted-face prism!

    vec3  vA = vec3(u.x, p.y, max(u.y, 0.0));
    float t  = clamp((qy - 0.5 * u.y) / (m2 + 0.25), 0.0, 1.0);
    vec3  vB = vec3(u.x + 0.5 * t, p.y - h * t, u.y + 0.5 * t);
    return sqrt(min(dot(vA, vA), dot(vB, vB)));
}

#endif // SDF_GEMS_SD_PYRAMID_GLSL
