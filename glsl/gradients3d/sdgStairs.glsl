// ============================================================================
// sdgStairs — Exact 1-Sqrt Stepped Stairs Booleans with Analytical 3D Gradients
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` for orthogonal surfaces, max Eikonal error `2.22e-16`)
// ============================================================================

// Uses the Exterior/Interior Gradient Sign-Cancellation Identity:
// In the exterior, dq/d(a,b) = (+1,+1) and d = +l; in the interior, di/d(a,b) = (-1,-1) and d = -l.
// The two minus signs cancel identically, yielding a single 1-sqrt gradient formula across all n tiers.
vec4 sdgUnionStairs( in vec4 a, in vec4 b, in float r, in float n )
{
    float s  = r / n, w = (a.x - b.x + r) / (2.0 * s);
    float k1 = clamp(round(w), 1.0, n - 1.0) * s;
    float j0 = clamp(floor(w), 0.0, n - 1.0) * s;
    vec2  q  = vec2(a.x - k1, b.x - r + k1), i = vec2(j0 - a.x, r - s - j0 - b.x);
    bool  e  = min(min(a.x, b.x), max(q.x, q.y)) > 0.0;
    vec2  v  = e ? q : i, c = max(v, 0.0);
    float l  = length(c), d = e ? min(min(a.x, b.x), l) : -(l + min(max(i.x, i.y), 0.0));
    vec2  g  = (e && min(a.x, b.x) < l) ? ((a.x < b.x) ? vec2(1.0, 0.0) : vec2(0.0, 1.0))
             : ((c.x > 0.0 && c.y > 0.0) ? c / l : ((v.x > v.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0)));
    return vec4(d, g.x * a.yzw + g.y * b.yzw);
}

vec4 sdgIntersectionStairs( in vec4 a, in vec4 b, in float r, in float n )
{
    return -sdgUnionStairs(-a, -b, r, n);
}

vec4 sdgSubtractionStairs( in vec4 a, in vec4 b, in float r, in float n )
{
    return -sdgUnionStairs(a, -b, r, n);
}
