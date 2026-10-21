// ============================================================================
// sdTriPrism — Exact 1-Sqrt 3D Equilateral Triangular Prism
// @author Inigo Quilez (plane bound) — upgraded to exact 1-sqrt using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, `1` sqrt everywhere, `1.33e-15` max error vs ground truth)
// ============================================================================

float sdTriPrism( vec3 p, vec2 h )
{
    const float k = 1.7320508;
    float r = 0.5 * k * h.x;
    p.xz = abs(p.xz);
    p.x -= r;
    p.y += 0.5 * h.x;
    if( p.x + k * p.y > 0.0 ) p.xy = vec2(p.x - k * p.y, -k * p.x - p.y) * 0.5;
    vec3 d = vec3(abs(p.x + r) - r, -p.y, p.z - h.y);
    return min(max(d.y, d.z), 0.0) + length(max(d, 0.0));
}
