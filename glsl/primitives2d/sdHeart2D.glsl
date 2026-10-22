// ============================================================================
// sdHeart2D & sdgHeart2D — Exact 0/1-Sqrt 2D Heart & Analytical Unit Gradient
// @author Inigo Quilez — upgraded to 0/1-sqrt using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in 24.4% of [-2,2]^2 for sdHeart2D, 0 NaNs on p.x == p.y)
// ============================================================================

// 1. 0/1-Sqrt Exact 2D Heart:
//    Uses the algebraic identity ||p - (0,1)||^2 - 0.5*(p.x - p.y)^2 == 0.5*(p.x + p.y - 1)^2 + 0.5 + (p.x - p.y).
//    Whenever 0 < p.x + p.y <= 1 and 0.5*(s - 1)^2 + 0.5 + d >= 0 (24.4% of [-2,2]^2),
//    the exact signed distance is simply (p.x - p.y) * sqrt(0.5) with 0 sqrts!
float sdHeart2D( in vec2 p )
{
    p.x = abs(p.x);
    float s = p.x + p.y, d = p.x - p.y;
    if( s > 1.0 ) return length(p - vec2(0.25, 0.75)) - 0.35355339;
    if( s > 0.0 && 0.5 * (s - 1.0) * (s - 1.0) + 0.5 + d >= 0.0 ) return d * 0.70710678;
    return (s <= 0.0) ? length(p) : -length(p - vec2(0.0, 1.0));
}

// 2. 0/1-Sqrt 2D Heart + Analytical Gradient (fixes 0.0 / 0.0 = NaN on p.x == p.y in IQ's sdgHeart)
vec3 sdgHeart2D( in vec2 p )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float s = p.x + p.y, d = p.x - p.y;
    if( s > 1.0 )
    {
        vec2  q = p - vec2(0.25, 0.75);
        float l = length(q);
        return vec3(l - 0.35355339, sx * q.x / l, q.y / l);
    }
    if( s > 0.0 && 0.5 * (s - 1.0) * (s - 1.0) + 0.5 + d >= 0.0 )
        return vec3(d * 0.70710678, sx * 0.70710678, -0.70710678);
    if( s <= 0.0 )
    {
        float l = max(length(p), 1e-12);
        return vec3(l, sx * p.x / l, p.y / l);
    }
    vec2  q = p - vec2(0.0, 1.0);
    float l = max(length(q), 1e-12);
    return vec3(-l, -sx * q.x / l, -q.y / l);
}
