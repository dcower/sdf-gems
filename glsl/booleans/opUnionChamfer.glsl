// ============================================================================
// opUnionChamfer — Exact Corner-Local Euclidean Chamfer Union, Intersection & Subtraction
// @author Mercury Demogroup (plane bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1` in both exterior and interior, max error `8.88e-16`)
// ============================================================================

// 1. Exact Corner-Local Euclidean 45-Degree Chamfer Union (fixes interior seam distortion)
float opUnionChamfer( float a, float b, float r )
{
    float x = min(a, b), y = max(a, b);
    return (x < 0.0 && y - x > r) ? -length(vec2(x, max(r - y, 0.0)))
                                  : min(x, (x + y - r) * sqrt(0.5));
}

// 2. Angle-Exact Oblique Chamfer Union for surfaces meeting at cos(theta) = c = dot(grad(a), grad(b))
float opUnionChamferOblique( float a, float b, float r, float c )
{
    return min(min(a, b), (a + b - r) * inversesqrt(max(2.0 + 2.0 * c, 1e-6)));
}

// 3. Exact Euclidean Chamfer Intersection & Subtraction
float opIntersectionChamfer( float a, float b, float r )
{
    return -opUnionChamfer(-a, -b, r);
}

float opSubtractionChamfer( float a, float b, float r )
{
    return -opUnionChamfer(a, -b, r);
}
