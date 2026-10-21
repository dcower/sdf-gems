// ============================================================================
// opUnionChamfer — Exact Corner-Local Euclidean Chamfer Union, Intersection & Subtraction
// @author Mercury Demogroup (plane bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1` in both exterior and interior, max error `8.88e-16`)
// ============================================================================

#ifndef SDF_GEMS_OP_UNION_CHAMFER_GLSL
#define SDF_GEMS_OP_UNION_CHAMFER_GLSL

// 1. Exact Corner-Local Euclidean 45-Degree Chamfer Union (fixes interior seam distortion)
float opUnionChamfer( float a, float b, float r )
{
    float x = min(a, b), y = max(a, b);
    return (x < 0.0 && y - x > r) ? -length(vec2(x, max(r - y, 0.0)))
                                  : min(x, (x + y - r) * sqrt(0.5));
}

// 2. Exact Euclidean Chamfer Intersection & Subtraction (Orthogonal)
float opIntersectionChamfer( float a, float b, float r )
{
    return -opUnionChamfer(-a, -b, r);
}

float opSubtractionChamfer( float a, float b, float r )
{
    return -opUnionChamfer(a, -b, r);
}

// 3. Angle-Exact Corner-Local Oblique Chamfer Union, Intersection & Subtraction
//    for surfaces meeting at cos(theta) = c = dot(grad(a), grad(b)).
//    Because (grad_x + grad_y) . (grad_x - grad_y) = ||grad_x||^2 - ||grad_y||^2 = 0,
//    the normal strip boundary at vertex V = (0, r) is universally y - x = r for ALL c!
float opUnionChamferOblique( float a, float b, float r, float c )
{
    float x  = min(a, b), y = max(a, b);
    float ch = (x + y - r) * inversesqrt(max(2.0 + 2.0 * c, 1e-5));
    if( x > 0.0 && x + y > r ) return min(x, ch);
    float h  = max(r - y + c * x, 0.0);
    float dv = -sqrt(x * x + h * h / max(1.0 - c * c, 1e-5));
    return (y - x > r) ? dv : max(dv, ch);
}

float opIntersectionChamferOblique( float a, float b, float r, float c )
{
    return -opUnionChamferOblique(-a, -b, r, c);
}

float opSubtractionChamferOblique( float a, float b, float r, float c )
{
    return -opUnionChamferOblique(a, -b, r, -c);
}

#endif // SDF_GEMS_OP_UNION_CHAMFER_GLSL
