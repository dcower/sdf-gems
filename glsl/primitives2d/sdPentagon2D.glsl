// ============================================================================
// sdPentagon2D, sdgPentagon2D & sdPentagonalPrism — 0/1-Sqrt 2D Pentagon, Gradient & 0/1-Sqrt 3D Prism
// @author Inigo Quilez (2D) — upgraded to 0/1-sqrt, 0/1-SFU gradient & 0/1-sqrt 3D prism using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & face strips; `0/1` sqrt in 3D prism)
// ============================================================================

#ifndef SDF_GEMS_SD_PENTAGON2D_GLSL
#define SDF_GEMS_SD_PENTAGON2D_GLSL

// 1. 0/1-Sqrt 2D Regular Pentagon (inradius r)
float sdPentagon2D( in vec2 p, in float r )
{
    const vec3 k = vec3(0.809016994, 0.587785252, 0.726542528);
    p.x = abs(p.x);
    p -= 2.0 * min(dot(vec2(-k.x, k.y), p), 0.0) * vec2(-k.x, k.y);
    p -= 2.0 * min(dot(vec2( k.x, k.y), p), 0.0) * vec2( k.x, k.y);
    vec2 d = vec2(abs(p.x) - r * k.z, p.y - r);
    return (d.x > 0.0) ? length(d) : d.y;
}

// 2. 0/1-SFU Analytical 2D Distance & Gradient vec3(d, grad.xy)
vec3 sdgPentagon2D( in vec2 p, in float r )
{
    const vec3 k  = vec3(0.809016994, 0.587785252, 0.726542528);
    const vec2 n1 = vec2(-0.809016994, 0.587785252);
    const vec2 n2 = vec2( 0.809016994, 0.587785252);
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float d1 = dot(n1, p); bool r1 = (d1 < 0.0); if( r1 ) p -= 2.0 * d1 * n1;
    float d2 = dot(n2, p); bool r2 = (d2 < 0.0); if( r2 ) p -= 2.0 * d2 * n2;

    vec2  v = vec2(abs(p.x) - r * k.z, p.y - r);
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
    if( r2 ) g -= 2.0 * dot(n2, g) * n2;
    if( r1 ) g -= 2.0 * dot(n1, g) * n1;
    g.x *= sx;
    return vec3(d, g);
}

// 3. Exact 0/1-Sqrt 3D Regular Pentagonal Prism (inradius h.x, half-height h.y)
float sdPentagonalPrism( in vec3 p, in vec2 h )
{
    const vec3 k = vec3(0.809016994, 0.587785252, 0.726542528);
    p.x = abs(p.x);
    p.xy -= 2.0 * min(dot(vec2(-k.x, k.y), p.xy), 0.0) * vec2(-k.x, k.y);
    p.xy -= 2.0 * min(dot(vec2( k.x, k.y), p.xy), 0.0) * vec2( k.x, k.y);
    vec3 d = vec3(abs(p.x) - h.x * k.z, p.y - h.x, abs(p.z) - h.y);
    return (max(d.x, min(d.y, d.z)) <= 0.0) ? max(d.y, d.z) : length(max(d, 0.0));
}

#endif // SDF_GEMS_SD_PENTAGON2D_GLSL
