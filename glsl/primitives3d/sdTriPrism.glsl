// ============================================================================
// sdTriPrism — 3D Triangular Prism (Plane Bound)
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  No (Plane Bound)
// ============================================================================

float sdTriPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);
    return max(q.z - h.y, max(q.x*0.866025 + p.y*0.5, -p.y) - h.x*0.5);
}
