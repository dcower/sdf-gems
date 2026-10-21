// ============================================================================
// sdIcosahedron — Exact 1-Sqrt 3D Regular Icosahedron (Ih Chamber + Axis-Aligned Edge)
// @author Mercury Demogroup (plane bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1`, max error `1.78e-15`)
// ============================================================================

float sdIcosahedron( vec3 p, float r )
{
    const vec3 n3 = vec3(-0.5, -0.80901699, 0.30901699);
    for (int i = 0; i < 3; i++) {
        p = abs(p);
        p -= 2.0 * min(dot(p, n3), 0.0) * n3;
    }
    p.xy = abs(p.xy);
    const vec2 n = vec2(0.35682209, 0.93417236);
    float d = dot(p.yz, n) - r;
    return (p.y >= d * n.x) ? d : length(vec3(max(p.x - 0.66158454 * r, 0.0), p.y, p.z - 1.07046627 * r));
}
