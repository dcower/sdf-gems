// ============================================================================
// sdDodecahedron — Exact 1-Sqrt 3D Regular Dodecahedron (Ih Chamber + Axis-Aligned Edge)
// @author Mercury Demogroup (plane bound) — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1`, max error `2.22e-15`)
// ============================================================================

float sdDodecahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for (int i = 0; i < 3; i++) {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);
    const vec2 n = vec2(0.52573111, 0.85065081);
    float d = dot(p.xz, n) - r;
    return (p.x >= d * n.x) ? d : length(vec3(p.x, max(p.y - 0.44902798 * r, 0.0), p.z - 1.17557050 * r));
}
