// ============================================================================
// sdTruncatedOctahedron — 3D Truncated Octahedron / Kelvin Cell (GDF Plane Bound)
// @author Mercury Demogroup
// @source https://mercury.sexy/hg_sdf/
// @exact  No (Plane Bound)
// ============================================================================

float sdTruncatedOctahedron( vec3 p, float r )
{
    p = abs(p);
    return max(max(p.x, max(p.y, p.z)) - r, (p.x + p.y + p.z - 1.5*r) * 0.57735027);
}
