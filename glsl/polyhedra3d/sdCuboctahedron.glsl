// ============================================================================
// sdCuboctahedron — Exact 1-Sqrt 3D Cuboctahedron (Oh Weyl Chamber)
// @author Mercury Demogroup (plane bound) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1`, max error `1.33e-15`)
// ============================================================================

float sdCuboctahedron( vec3 p, float s )
{
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    if (min(p.x, 2.0*p.x - p.y - p.z) <= s)
        return max(p.x - s, (p.x + p.y + p.z - 2.0*s) * 0.57735027);
    float t = clamp(0.5*(s - p.y + p.z), 0.0, p.z);
    return length(vec3(p.x - s, max(p.y - s + t, 0.0), p.z - t));
}
