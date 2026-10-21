// ============================================================================
// opRepBCC — 1-Tap 3D Body-Centered Cubic (BCC) Kelvin-Cell Repetition
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (symmetric SDFs)
// ============================================================================

vec3 opRepBCC( in vec3 p, in float s )
{
    vec3 q  = p / s;
    vec3 p0 = round(q), r = q - p0;
    vec3 c  = (abs(r.x) + abs(r.y) + abs(r.z) > 0.75) ? (round(q - 0.5) + 0.5) : p0;
    return (q - c) * s;
}
