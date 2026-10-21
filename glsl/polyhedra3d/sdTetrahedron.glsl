// ============================================================================
// sdTetrahedron — 3D Regular Tetrahedron (GDF Plane Bound)
// @author Mercury Demogroup
// @source https://mercury.sexy/hg_sdf/
// @exact  No (Plane Bound)
// ============================================================================

float sdTetrahedron( vec3 p, float r )
{
    return (max(abs(p.x + p.y) - p.z, abs(p.x - p.y) + p.z) - r) * 0.57735027;
}
