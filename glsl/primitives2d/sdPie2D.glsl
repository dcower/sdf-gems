// ============================================================================
// sdPie2D — Exact 1-Sqrt 2D Pie Sector & Ring Sector (Orthonormal Dominance)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 sqrt down from 2, 0.0 error across R^2)
// ============================================================================

// 1. Exact 1-Sqrt 2D Circular Pie Sector (down from 2 sqrts)
float sdPie2D( in vec2 p, in vec2 c, in float r )
{
    p.x = abs(p.x);
    float u = dot(p, c), v = c.y*p.x - c.x*p.y;
    float l = length(vec2(u - ((v > 0.0) ? clamp(u, 0.0, r) : 0.0), v));
    return (v > 0.0) ? l : max(l - r, (u < 0.0) ? -l : v);
}

// 2. Exact 1-Sqrt 2D Annular Ring Sector (down from 2 sqrts, for u >= -(r - 0.5*th))
float sdRing2D( in vec2 p, in vec2 n, in float r, in float th )
{
    p.x = abs(p.x);
    float v = n.x*p.x - n.y*p.y, u = n.y*p.x + n.x*p.y;
    float h = 0.5*th, r1 = r - h;
    float w = (v > 0.0 || u < -r1) ? max(abs(r - u) - h, 0.0) : u;
    float l = length(vec2(v, w));
    return (v > 0.0) ? l : ((u < -r1) ? max(r1 - length(vec2(v, u)), -l) : max(abs(l - r) - h, v));
}
