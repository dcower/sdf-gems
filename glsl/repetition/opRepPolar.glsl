// ============================================================================
// opRepPolar — Gap-Aware 2-Tap Limited Polar Repetition (opLimitedPolarRep2Tap)
// @author Mercury Demogroup — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (eliminates theta = +-pi branch cut)
// ============================================================================

void opLimitedPolarRep2Tap( in vec2 p, in float n, in float max_idx, out float id0, out float id1 )
{
    float sp = SDF_TAU / n;
    float a  = atan(p.y, p.x) / sp;
    float i0 = clamp(floor(a), -max_idx, max_idx);
    float i1 = clamp(i0 + 1.0, -max_idx, max_idx);
    if (abs(a) > max_idx) { i0 = -max_idx; i1 = max_idx; }
    id0 = i0 * sp;
    id1 = i1 * sp;
}
