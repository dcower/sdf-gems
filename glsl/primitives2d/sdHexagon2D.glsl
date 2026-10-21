// ============================================================================
// sdHexagon2D — Exact 0/1-Sqrt 2D Regular Hexagon (Fundamental Chamber Vertex Dominance)
// @author Inigo Quilez — upgraded to 0-sqrt interior/face-strip using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in interior & face strips, `1` in corner wedges)
// ============================================================================

#ifndef SDF_GEMS_SD_HEXAGON2D_GLSL
#define SDF_GEMS_SD_HEXAGON2D_GLSL

float sdHexagon2D( vec2 p, float r )
{
    const vec3 k = vec3(-0.866025404, 0.5, 0.577350269);
    p = abs(p);
    p -= 2.0*min(dot(k.xy, p), 0.0)*k.xy;
    vec2 d = vec2(abs(p.x) - k.z*r, p.y - r);
    return (d.x > 0.0) ? length(d) : d.y;
}

#endif // SDF_GEMS_SD_HEXAGON2D_GLSL
