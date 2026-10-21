// ============================================================================
// sdOctahedron — Exact 3D Regular Octahedron (4-Line Branchless Pythagorean Weyl-Basis Form)
// @author Inigo Quilez — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 0 branches)
// ============================================================================

float sdOctahedron( vec3 p, float s )
{
    p = abs(p);
    float m = p.x + p.y + p.z - s, z = min(p.x, min(p.y, p.z));
    float a = max(m - 3.0*z, 0.0);
    float b = max(2.0*max(p.x, max(p.y, p.z)) + z - m - 2.0*s, 0.0);
    return sign(m) * sqrt((2.0*m*m + a*a + 3.0*b*b) / 6.0);
}
