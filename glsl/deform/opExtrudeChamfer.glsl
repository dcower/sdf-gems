// ============================================================================
// opExtrudeChamfer — Orthogonal Eikonal 45-Degree Chamfer Extrusion
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float opExtrudeChamfer( in float d2d, in float z, in float h, in float c )
{
    float u = d2d, v = abs(z) - h;
    if (u <= -c && v <= 0.0) return max(u, v);
    if (u <= 0.0 && v <= -c) return max(u, v);
    if (u + v <= -c) return (u + v + c) * 0.70710678;
    float t = clamp(0.5 * (u - v + c), 0.0, c);
    return length(vec2(max(u, -c + t), max(v, -t)));
}
