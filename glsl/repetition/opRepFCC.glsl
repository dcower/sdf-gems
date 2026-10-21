// ============================================================================
// opRepFCC — 1-Tap 3D Face-Centered Cubic (FCC / D3) Rhombic-Dodecahedron Repetition
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (symmetric SDFs)
// ============================================================================

vec3 opRepFCC( in vec3 p, in float s )
{
    vec3  q = p / s;
    vec3  r = round(q);
    float m = mod(r.x + r.y + r.z, 2.0);
    if (m != 0.0) {
        vec3 d = abs(q - r);
        vec3 sgn = vec3(greaterThanEqual(q, r)) * 2.0 - 1.0;
        r += (d.x >= d.y && d.x >= d.z) ? vec3(sgn.x, 0.0, 0.0) :
             (d.y >= d.z)               ? vec3(0.0, sgn.y, 0.0) : vec3(0.0, 0.0, sgn.z);
    }
    return (q - r) * s;
}
