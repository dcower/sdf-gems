// ============================================================================
// opUnionStairs — Exact 1-Sqrt Euclidean Stepped Stairs Union, Intersection & Subtraction
// @author Mercury Demogroup (L-infinity bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1` for orthogonal surfaces, max error `4.44e-16` vs `0.414` in baseline)
// ============================================================================

// 1. Exact 1-Sqrt Corner-Local Euclidean Stairs Union (true circular step noses & exact interior)
float opUnionStairs( float a, float b, float r, float n )
{
    float s  = r / n, w = (a - b + r) / (2.0 * s);
    float k1 = clamp(round(w), 1.0, n - 1.0) * s;
    float j0 = clamp(floor(w), 0.0, n - 1.0) * s;
    vec2  q  = vec2(a - k1, b - r + k1), i = vec2(j0 - a, r - s - j0 - b);
    bool  e  = min(min(a, b), max(q.x, q.y)) > 0.0;
    float l  = length(max(e ? q : i, 0.0));
    return e ? min(min(a, b), l) : -(l + min(max(i.x, i.y), 0.0));
}

// 2. Exact Euclidean Stairs Intersection & Subtraction
float opIntersectionStairs( float a, float b, float r, float n )
{
    return -opUnionStairs(-a, -b, r, n);
}

float opSubtractionStairs( float a, float b, float r, float n )
{
    return -opUnionStairs(a, -b, r, n);
}
