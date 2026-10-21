// ============================================================================
// sdCone — Exact 1-to-2 Sqrt 3D Finite & Infinite Cone
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt` everywhere except base-rim wedge & 1 div, down from 2 `sqrt`s & 3 divs)
// ============================================================================

// c = vec2(sin(alpha), cos(alpha)) normalized cone aperture, h = cone height
float sdCone( vec3 p, vec2 c, float h )
{
    float x2 = dot(p.xz, p.xz);
    if( p.y >= 0.0 && c.x * c.x * x2 <= c.y * c.y * p.y * p.y )
        return sqrt(x2 + p.y * p.y);
    float wx = sqrt(x2), rx = wx - h * c.x / c.y, d = -p.y - h;
    return (rx > 0.0 && c.x * rx + c.y * d >= 0.0)
         ? length(vec2(rx, d))
         : max(c.y * wx + c.x * p.y, d);
}

// Infinite 3D cone (1 `sqrt` everywhere in R^3, down from 2 `sqrt`s)
float sdConeInf( vec3 p, vec2 c )
{
    float x2 = dot(p.xz, p.xz);
    return (p.y >= 0.0 && c.x * c.x * x2 <= c.y * c.y * p.y * p.y)
         ? sqrt(x2 + p.y * p.y)
         : sqrt(x2) * c.y + p.y * c.x;
}
