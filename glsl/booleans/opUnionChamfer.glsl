// ============================================================================
// opUnionChamfer — 45-Degree Chamfer Union
// @author Mercury Demogroup
// @source https://mercury.sexy/hg_sdf/
// @exact  Bound
// ============================================================================

float opUnionChamfer( float a, float b, float r )
{
    return min(min(a, b), (a - r + b)*sqrt(0.5));
}
