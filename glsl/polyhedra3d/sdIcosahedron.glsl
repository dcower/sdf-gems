// ============================================================================
// sdIcosahedron — 3D Regular Icosahedron (GDF Plane Bound)
// @author Mercury Demogroup
// @source https://mercury.sexy/hg_sdf/
// @exact  No (Plane Bound)
// ============================================================================

float sdIcosahedron( vec3 p, float r )
{
    p = abs(p);
    const vec2 n = vec2(0.35682209, 0.93417236);
    float d1 = (p.x + p.y + p.z) * 0.57735027;
    float d2 = max(dot(p.xy, n), max(dot(p.yz, n), dot(p.zx, n)));
    return max(d1, d2) - r;
}
