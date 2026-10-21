// ============================================================================
// opUnionStairs — Stepped Architectural Stairs Union
// @author Mercury Demogroup
// @source https://mercury.sexy/hg_sdf/
// @exact  Bound
// ============================================================================

float opUnionStairs( float a, float b, float r, float n )
{
    float s = r/n;
    float u = b - r;
    return min(min(a, b), 0.5*(u + a + abs((mod(u - a + s, 2.0*s)) - s)));
}
