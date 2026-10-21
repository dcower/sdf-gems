// ============================================================================
// opExtrudeRound — Orthogonal Eikonal Round-Bevel Extrusion
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float opExtrudeRound( in float d2d, in float z, in float h, in float r )
{
    vec2 q = vec2(d2d + r, abs(z) - h + r);
    return min(max(q.x, q.y), 0.0) + length(max(q, 0.0)) - r;
}
