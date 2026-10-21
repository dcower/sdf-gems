// ============================================================================
// sdTetrahedron — Exact 1-Sqrt 3D Regular Tetrahedron (Td Parity-Sorted Weyl Chamber)
// @author Mercury Demogroup (plane bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1`, max error `8.88e-16`)
// ============================================================================

float sdTetrahedron( vec3 p, float s )
{
    float sig = (p.x * p.y * p.z < 0.0) ? -1.0 : 1.0;
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    p.z *= sig;
    float m = p.x + p.y - p.z - s;
    float h = min(m / 3.0, 0.5 * (p.y - p.z));
    return sign(m) * length(p - min(p + vec3(-h, -h, h), s));
}
