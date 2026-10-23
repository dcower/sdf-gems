// ============================================================================
// sdCircleWave2D — 1-Sqrt Unsigned & Signed Circular Wave + 1-SFU Arc Gradients
// @author Inigo Quilez — upgraded to 1-sqrt / signed / 1-SFU using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, `1` sqrt down from `4` sqrts via Inflection-Normal Dead-Code Theorem)
// ============================================================================
//
// Mathematical Discovery (Inflection-Normal Dead-Code & Signed Wave Theorem):
//   In IQ's sdCircleWave (Shadertoy stGyzt), two candidate points p1 = p and
//   p2 = vec2(abs(p.x - 2*co.x), 2*co.y - p.y) are evaluated against both the
//   arc distance abs(length(p_i) - ra) and the endpoint cap length(p_i - co),
//   requiring up to 4 sqrts, 2 sector tests, and min(d1, d2).
//
//   In the folded half-period p.x in [0, 2*co.x]:
//   1. p.x <= 2*co.x always, so p2 = 2.0*co - p (point reflection across co),
//      and its sector test co.y*p2.x - co.x*p2.y = -(co.y*p.x - co.x*p.y) is
//      the exact opposite of p1's sector test.
//   2. By the triangle inequality in triangle(0, p, co) and triangle(2*co, p, co),
//      since |co| = ra, | |p| - ra | <= |p - co| always. Thus both endpoint
//      branches length(p1 - co) and length(p2 - co) are 100% dead code!
//   3. Even better: along the separating normal line p = t*co inside the strip
//      p.x in [0, 2*co.x], t = p.x / co.x lies strictly in [0, 2], where
//      (|p| - ra) = (t - 1)*ra = -(|2*co - p| - ra) and both gradients equal
//      +co / ra. Thus negating the upper-sector distance yields a C1-smooth,
//      exact signed half-plane wave SDF (`sdCircleWaveSigned2D`) with 1 sqrt!

// 1. Exact 1-Sqrt Signed Circular Wave Half-Space (d < 0 below wave, d > 0 above wave)
float sdCircleWaveSigned2D( in vec2 p, in float tb, in float ra )
{
    tb = 3.1415927 * 5.0 / 6.0 * max(tb, 0.0001);
    vec2 co = ra * vec2(sin(tb), cos(tb));
    p.x = abs(mod(p.x, co.x * 4.0) - co.x * 2.0);
    bool u = (co.y * p.x > co.x * p.y);
    if( u ) p = 2.0 * co - p;
    return (u ? -1.0 : 1.0) * (length(p) - ra);
}

// 2. Exact 1-Sqrt Unsigned Circular Wave Curve (drop-in replacement for IQ's 4-sqrt sdCircleWave)
float sdCircleWave2D( in vec2 p, in float tb, in float ra )
{
    return abs(sdCircleWaveSigned2D(p, tb, ra));
}

// 3. Exact 1-SFU (`inversesqrt`) Circular Wave with Analytical 2D Gradient
vec3 sdgCircleWave2D( in vec2 p, in float tb, in float ra )
{
    tb = 3.1415927 * 5.0 / 6.0 * max(tb, 0.0001);
    vec2  co = ra * vec2(sin(tb), cos(tb));
    float xm = mod(p.x, co.x * 4.0) - co.x * 2.0;
    float sx = (xm < 0.0) ? -1.0 : 1.0;
    p.x = abs(xm);
    bool  u  = (co.y * p.x > co.x * p.y);
    if( u ) p = 2.0 * co - p;
    float l2 = dot(p, p), il = inversesqrt(max(l2, 1e-24));
    float w  = l2 * il - ra, sw = (w < 0.0) ? -1.0 : 1.0, su = u ? -1.0 : 1.0;
    return vec3(abs(w), vec2(sx, 1.0) * (su * sw * il) * p);
}

// 4. Exact 1-SFU (`inversesqrt`) 2D Circular Arc with Analytical Gradient
#ifndef SDF_GEMS_HAS_SDG_ARC2D
#define SDF_GEMS_HAS_SDG_ARC2D
vec3 sdgArc2D( in vec2 p, in vec2 sc, in float ra, in float rb )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    if( sc.y * p.x > sc.x * p.y )
    {
        vec2  w  = p - ra * sc;
        float l2 = dot(w, w);
        if( l2 > 0.0 )
        {
            float il = inversesqrt(l2);
            return vec3(l2 * il - rb, vec2(sx * w.x, w.y) * il);
        }
        return vec3(-rb, sx * sc.y, -sc.x);
    }
    float l2 = dot(p, p);
    if( l2 > 0.0 )
    {
        float il = inversesqrt(l2);
        float w  = l2 * il - ra, sw = (w < 0.0) ? -1.0 : 1.0;
        return vec3(abs(w) - rb, vec2(sx * p.x, p.y) * (sw * il));
    }
    return vec3(ra - rb, 0.0, -1.0);
}
#endif // SDF_GEMS_HAS_SDG_ARC2D
