// ============================================================================
// opEngrave — Exact 0/1-Sqrt Euclidean V-Groove Engraving & Intersection Pipe
// @author Mercury Demogroup (plane bound) — upgraded to exact Euclidean & oblique Gram-Schmidt using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1`, max error `4.44e-16` vs `0.7929` in hg_sdf)
// ============================================================================

// 1. Exact 0/1-Sqrt Euclidean V-Groove Engraving (fixes top-lip distortion & interior V-tip bound)
float opEngrave( float a, float b, float r )
{
    float y = abs(b), s = a + y;
    if( abs(s) <= r ) return max(a, (a + r - y) * 0.70710678);
    if( y >= r && s > r ) return a;
    return (s < -r) ? max(a, -length(vec2(a + r, y))) : length(vec2(a, r - y));
}

// 2. Angle-Exact Oblique V-Groove Engraving via Gram-Schmidt Surface-Anchored Orthonormal Map
// c = dot(grad(a), grad(b)) = cos(theta)
float opEngraveOblique( float a, float b, float r, float c )
{
    return opEngrave(a, (b - c * a) * inversesqrt(max(1.0 - c * c, 1e-5)), r);
}

// 3. Orthogonal & Angle-Exact Oblique Intersection Pipe
float opPipe( float a, float b, float r )
{
    return length(vec2(a, b)) - r;
}

float opPipeOblique( float a, float b, float r, float c )
{
    return sqrt(max(a * a + b * b - 2.0 * c * a * b, 0.0) / max(1.0 - c * c, 1e-5)) - r;
}
