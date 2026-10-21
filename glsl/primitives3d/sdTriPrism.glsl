// ============================================================================
// sdTriPrism — Exact 1-Sqrt 3D Equilateral Triangular Prism
// @author Inigo Quilez (plane bound) — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, upgraded from plane bound)
// ============================================================================

float sdTriPrism( vec3 p, vec2 h )
{
    const vec2 k = vec2(0.5, 0.8660254);
    p.xz = abs(p.xz);
    p.xy -= 2.0*max(dot(k, p.xy), 0.0)*k;
    vec3 d = vec3(abs(p.x) - k.y*h.x, -p.y - 0.5*h.x, p.z - h.y);
    return min(max(d.y, d.z), 0.0) + length(max(d, 0.0));
}
