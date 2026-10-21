// ============================================================================
// opExtrudeStairs — Exact 1-Sqrt Orthogonal Eikonal n-Tier Ziggurat Extrusion
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @exact  Yes (`||grad d|| = 1` across exterior circular step noses, treads, and interior notches)
// ============================================================================

float opExtrudeStairs( float d2d, float z, float h, float r, float n )
{
    float u  = d2d, v = abs(z) - h;
    float s  = r / n, w = (v - u + r) / (2.0 * s);
    float j0 = clamp(floor(w), 0.0, n - 1.0) * s;
    float k1 = clamp(round(w), 1.0, n - 1.0) * s;
    vec2  p  = vec2(u + j0, v + r - s - j0);
    vec2  q  = vec2(-u - k1, -v - r + k1);
    bool  m  = max(p.x, p.y) < 0.0;
    float l  = length(max(m ? q : p, 0.0));
    return m ? max(max(u, v), -l) : l;
}
