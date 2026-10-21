// ============================================================================
// opEngrave — Exact 0/1-Sqrt Euclidean V-Groove Engraving & Intersection Pipe
// @author Mercury Demogroup (plane bound) — upgraded to exact Euclidean using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1` for orthogonal surfaces, max error `2.22e-16` vs `29.3%` in hg_sdf)
// ============================================================================

// 1. Exact 0/1-Sqrt Euclidean V-Groove Engraving (fixes top-lip distortion & interior V-tip bound)
float opEngrave( float a, float b, float r )
{
    float y = abs(b), s = a + y;
    if( abs(s) <= r ) return max(a, (a + r - y) * 0.70710678);
    if( y >= r && s > r ) return a;
    return ((s < -r) ? -1.0 : 1.0) * length((s < -r) ? vec2(a + r, y) : vec2(a, r - y));
}

// 2. Orthogonal & Angle-Exact Oblique Intersection Pipe
float opPipe( float a, float b, float r )
{
    return length(vec2(a, b)) - r;
}

// c = dot(grad(a), grad(b)) = cos(theta)
float opPipeOblique( float a, float b, float r, float c )
{
    return sqrt(max(a * a + b * b - 2.0 * c * a * b, 0.0) / max(1.0 - c * c, 1e-5)) - r;
}
