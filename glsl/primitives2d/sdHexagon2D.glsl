// ============================================================================
// sdHexagon2D & sdgHexagon2D — Exact 0/1-Sqrt 2D Regular Hexagon & Analytical Gradient
// @author Inigo Quilez — upgraded to 0-sqrt interior/face-strip & 0/1-SFU gradient using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & face strips, `1` in corner wedges)
// ============================================================================

#ifndef SDF_GEMS_SD_HEXAGON2D_GLSL
#define SDF_GEMS_SD_HEXAGON2D_GLSL

float sdHexagon2D( vec2 p, float r )
{
    const vec3 k = vec3(-0.866025404, 0.5, 0.577350269);
    p = abs(p);
    p -= 2.0 * min(dot(k.xy, p), 0.0) * k.xy;
    vec2 d = vec2(abs(p.x) - k.z * r, p.y - r);
    return (d.x > 0.0) ? length(d) : d.y;
}

vec3 sdgHexagon2D( vec2 p, float r )
{
    const vec3 k = vec3(-0.866025404, 0.5, 0.577350269);
    vec2 s = vec2((p.x < 0.0) ? -1.0 : 1.0, (p.y < 0.0) ? -1.0 : 1.0);
    p = abs(p);
    float dk = dot(k.xy, p);
    bool  rf = (dk < 0.0);
    if( rf ) p -= 2.0 * dk * k.xy;

    vec2  v = vec2(abs(p.x) - k.z * r, p.y - r);
    float d;
    vec2  g;
    if( v.x <= 0.0 )
    {
        d = v.y;
        g = vec2(0.0, 1.0);
    }
    else
    {
        float v2   = dot(v, v);
        float invV = inversesqrt(max(v2, 1e-30));
        d = v2 * invV;
        g = vec2((p.x < 0.0) ? -v.x : v.x, v.y) * invV;
    }
    if( rf ) g -= 2.0 * dot(k.xy, g) * k.xy;
    return vec3(d, s * g);
}

#endif // SDF_GEMS_SD_HEXAGON2D_GLSL
