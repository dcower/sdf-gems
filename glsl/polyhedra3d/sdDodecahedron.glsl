// ============================================================================
// sdDodecahedron — 3D Regular Dodecahedron (GDF Plane Bound)
// @author Mercury Demogroup
// @source https://mercury.sexy/hg_sdf/
// @exact  No (Plane Bound)
// ============================================================================

float sdDodecahedron( vec3 p, float r )
{
    p = abs(p);
    const vec2 n = vec2(0.52573111, 0.85065081);
    return max(dot(p.xy, n), max(dot(p.yz, n), dot(p.zx, n))) - r;
}
