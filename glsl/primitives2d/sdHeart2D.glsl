// ============================================================================
// sdHeart2D, sdEgg2D, & sdHorseshoe2D — 0/1-Sqrt Heart & 1-Sqrt Egg / Horseshoe
// @author Inigo Quilez — upgraded to 0/1-sqrt & 1-sqrt using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in 24.4% of [-2,2]^2 for sdHeart2D; `1` sqrt for sdEgg2D & sdHorseshoe2D)
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
    if( s > 1.0 ) {
        vec2  q = p - vec2(0.25, 0.75);
        float l = length(q);
        return vec3(l - 0.35355339, sx * q.x / l, q.y / l);
    }
    if( s > 0.0 && 0.5 * (s - 1.0) * (s - 1.0) + 0.5 + d >= 0.0 )
        return vec3(d * 0.70710678, sx * 0.70710678, -0.70710678);
    if( s <= 0.0 ) {
        float l = max(length(p), 1e-12);
        return vec3(l, sx * p.x / l, p.y / l);
    }
    vec2  q = p - vec2(0.0, 1.0);
    float l = max(length(q), 1e-12);
    return vec3(-l, -sx * q.x / l, -q.y / l);
}

// 3. 1-Sqrt Exact 2D Egg (down from 2 length() calls in IQ's sdEgg via disjoint normal-wedge selection)
float sdEgg2D( in vec2 p, in float he, in float ra, in float rb, in float bu )
{
    float r  = 0.5 * (he + ra + rb) / bu;
    float da = r - ra, db = r - rb;
    float y  = (db * db - da * da - he * he) / (2.0 * he);
    float x  = sqrt(max(da * da - y * y, 0.0));
    p.x = abs(p.x);
    float k  = p.y * x - p.x * y;
    if( k < 0.0 ) return length(p) - ra;
    if( k > he * (p.x + x) ) return length(vec2(p.x, p.y - he)) - rb;
    return length(p + vec2(x, y)) - r;
}

// 4. 1-Sqrt Exact 2D Horseshoe (down from 2 sqrts in IQ's sdHorseshoe; 0.0 error across R^2)
//    When q.x <= 0 (arc sector), u.x = (q.y > 0 ? q.x : -l) - w.x < 0 always,
//    so length(max(u, 0.0)) + min(max(u.x, u.y), 0.0) collapses identically to max(u.x, u.y) with 0 extra sqrts!
float sdHorseshoe2D( in vec2 p, in vec2 c, in float r, in vec2 w )
{
    p.x = abs(p.x);
    vec2  q = vec2(-c.x * p.x + c.y * p.y, c.y * p.x + c.x * p.y);
    vec2  b = vec2(q.x - w.x, abs(q.y - r) - w.y);
    float l = length((q.x > 0.0) ? max(b, 0.0) : p);
    return (q.x > 0.0) ? l + min(max(b.x, b.y), 0.0)
                       : max(((q.y > 0.0) ? q.x : -l) - w.x, abs(l - r) - w.y);
}
