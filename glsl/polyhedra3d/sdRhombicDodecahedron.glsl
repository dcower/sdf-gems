// ============================================================================
// sdRhombicDodecahedron — Exact 4-Line Branchless 3D Rhombic Dodecahedron (FCC Voronoi Cell)
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, max error `8.88e-16`)
// ============================================================================

float sdRhombicDodecahedron( vec3 p, float s )
{
    p = abs(p);
    float x = max(p.x, max(p.y, p.z)), z = min(p.x, min(p.y, p.z));
    float y = p.x + p.y + p.z - x - z, m = x + y - s;
    float a = max(x - y + 2.0*z - s, 0.0);
    float b = max(abs(y + z - x + 0.25*s) - 0.75*s, 0.0);
    return sign(m) * sqrt((3.0*m*m + a*a + 2.0*b*b) / 6.0);
}
