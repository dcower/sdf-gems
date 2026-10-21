// ============================================================================
// opRepPolar — Polar / Rotational Domain Repetition
// @author Mercury Demogroup
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (symmetric SDFs)
// ============================================================================

vec2 opRepPolar( in vec2 p, in float repetitions )
{
    float angle = SDF_TAU / repetitions;
    float a = atan(p.y, p.x) + angle*0.5;
    float r = length(p);
    a = mod(a, angle) - angle*0.5;
    return vec2(cos(a), sin(a))*r;
}
