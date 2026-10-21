// ============================================================================
// sdTruncatedOctahedron — Exact 1-Sqrt 3D Truncated Octahedron / Kelvin Cell
// @author Mercury Demogroup (plane bound) — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://mercury.sexy/hg_sdf/
// @exact  Yes (`||grad d|| = 1`, max error `1.78e-15`)
// ============================================================================

float sdTruncatedOctahedron( vec3 p, float s )
{
    p = abs(p);
    if (p.x < p.y) p.xy = p.yx;
    if (p.y < p.z) p.yz = p.zy;
    if (p.x < p.y) p.xy = p.yx;
    float d_hex = (p.x + p.y + p.z - 1.5*s) * 0.57735027;
    if ((p.x <= s && d_hex <= 0.0) || max(2.0*p.x - p.y - p.z, p.x + p.y - 2.0*p.z) <= 1.5*s)
        return max(p.x - s, d_hex);
    float t1 = clamp(0.5*(0.5*s - p.y + p.z), 0.0, p.z);
    vec3  d1 = vec3(p.x - s, max(p.y - 0.5*s + t1, 0.0), p.z - t1);
    float t2 = max(0.5*(0.5*s - p.x + p.y), 0.0);
    vec3  d2 = vec3(p.x - s + t2, p.y - 0.5*s - t2, p.z);
    return sqrt(min(dot(d1, d1), dot(d2, d2)));
}
