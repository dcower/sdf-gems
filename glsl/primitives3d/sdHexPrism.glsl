// ============================================================================
// sdHexPrism — Exact 3D Hexagonal Prism (1-Sqrt Fundamental Sector Optimization)
// @author Inigo Quilez — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdHexPrism( vec3 p, vec2 h )
{
    const vec3 k = vec3(-0.8660254, 0.5, 0.57735027);
    p = abs(p);
    p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
    vec3 d = vec3(abs(p.x) - k.z*h.x, p.y - h.x, p.z - h.y);
    return min(max(d.y, d.z), 0.0) + length(max(d, 0.0));
}
