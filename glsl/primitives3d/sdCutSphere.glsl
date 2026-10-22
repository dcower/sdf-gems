// ============================================================================
// sdCutSphere & sdgCutSphere — Exact 0/1/2-SFU 3D Cut Sphere & Cancellation-Free Rim Gradient
// @author Inigo Quilez — upgraded to unified 0-sqrt cap-Voronoi & 2-SFU rim gradient using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `0` SFUs in exterior + interior cap zone, `1` SFU in sphere, `2` SFUs in rim)
// ============================================================================

#ifndef SDF_GEMS_SD_CUT_SPHERE_GLSL
#define SDF_GEMS_SD_CUT_SPHERE_GLSL

float sdCutSphere( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), dy = h - p.y;
    if( x2 <= (r + h) * (r - h + 2.0 * min(dy, 0.0)) ) return dy;
    float w2 = r * r - h * h;
    if( h * abs(h) * x2 < w2 * p.y * abs(p.y) ) return sqrt(x2 + p.y * p.y) - r;
    float d = x2 - w2;
    return sqrt(d * d / (x2 + w2 + 2.0 * sqrt(w2 * x2)) + dy * dy);
}

vec4 sdgCutSphere( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), dy = h - p.y;
    if( x2 <= (r + h) * (r - h + 2.0 * min(dy, 0.0)) )
        return vec4(dy, 0.0, -1.0, 0.0);
    float w2 = r * r - h * h;
    if( h * abs(h) * x2 < w2 * p.y * abs(p.y) )
    {
        float r2 = x2 + p.y * p.y, invR = inversesqrt(max(r2, 1e-30));
        return vec4(r2 * invR - r, p * invR);
    }
    float d  = x2 - w2, wx = sqrt(w2 * x2);
    float sc = d / max(x2 + wx, 1e-30);
    vec3  v  = vec3(sc * p.x, -dy, sc * p.z);
    float v2 = dot(v, v), invV = inversesqrt(max(v2, 1e-30));
    return vec4(v2 * invV, v * invV);
}

#endif // SDF_GEMS_SD_CUT_SPHERE_GLSL
