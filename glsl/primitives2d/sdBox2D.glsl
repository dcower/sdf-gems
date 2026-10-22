// ============================================================================
// sdBox2D, sdRoundBox2D, sdgBox2D & sdgRoundBox2D — Exact 0/1-SFU 2D Box + Gradients
// @author Inigo Quilez — upgraded to 0-SFU interior & face-slab fast path using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 0 SFUs in interior + all 4 exterior face slabs, 1 in corners)
// ============================================================================

#ifndef SDF_GEMS_SD_BOX2D_GLSL
#define SDF_GEMS_SD_BOX2D_GLSL

float sdBox2D( vec2 p, vec2 b )
{
    vec2 d = abs(p) - b;
    return (min(d.x, d.y) <= 0.0) ? max(d.x, d.y) : length(d);
}

float sdRoundBox2D( vec2 p, vec2 b, float r )
{
    vec2 d = abs(p) - b + r;
    return ((min(d.x, d.y) <= 0.0) ? max(d.x, d.y) : length(d)) - r;
}

vec3 sdgBox2D( vec2 p, vec2 b )
{
    vec2 w = abs(p) - b;
    vec2 s = vec2(p.x < 0.0 ? -1.0 : 1.0, p.y < 0.0 ? -1.0 : 1.0);
    if( min(w.x, w.y) <= 0.0 )
        return (w.x > w.y) ? vec3(w.x, s.x, 0.0) : vec3(w.y, 0.0, s.y);
    float l2 = dot(w, w), invL = inversesqrt(l2);
    return vec3( l2 * invL, s * w * invL );
}

vec3 sdgRoundBox2D( vec2 p, vec2 b, float r )
{
    vec3 res = sdgBox2D( p, b - r );
    return vec3( res.x - r, res.yz );
}

#endif // SDF_GEMS_SD_BOX2D_GLSL
