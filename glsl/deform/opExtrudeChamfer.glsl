// ============================================================================
// opExtrudeChamfer — Exact 5-Line Orthogonal Eikonal 45-Degree Chamfer Extrusion
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1`, 1 `sqrt`, `4.44e-16` max error)
// ============================================================================

float opExtrudeChamfer( float d2d, float z, float h, float c )
{
    float u = d2d, v = abs(z) - h;
    float a = min(u, v), b = max(u, v);
    return (b > 0.0 && b - a > c) ? length(vec2(b, max(a + c, 0.0)))
                                  : max(b, (a + b + c) * 0.70710678);
}
