// ============================================================================
// sdCappedCylinder — 0/1/2-Sqrt Exact 3D Vertical Capped Cylinder
// @author Inigo Quilez — upgraded to 0/1/2 sqrts via Cap-Voronoi Radius Theorem using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `0` sqrts in cap zones, `1` in wall zone, `2` in rim wedge)
// ============================================================================

float sdCappedCylinder( vec3 p, float h, float r )
{
    float r2 = dot(p.xz, p.xz), wy = abs(p.y) - h, rc = r + min(wy, 0.0);
    if( rc >= 0.0 && r2 <= rc * rc ) return wy;
    float wx = sqrt(r2) - r;
    return (wy <= 0.0) ? wx : sqrt(wx * wx + wy * wy);
}
