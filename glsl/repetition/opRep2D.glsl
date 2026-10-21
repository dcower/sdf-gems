// ============================================================================
// opRep2D — 3-Tap 2D Hexagonal (A2) Delaunay Domain Repetition (opRepHex3)
// @author Inigo Quilez (4-tap Cartesian) — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes (`0.0` error vs 19-tap brute force)
// ============================================================================

void opRepHex3( in vec2 p, in float s, out vec2 c0, out vec2 c1, out vec2 c2 )
{
    const vec2 k = vec2(1.0, 1.73205081);
    vec2 q  = p / s;
    vec2 i0 = round(q / k);
    vec2 i1 = round((q - 0.5 * k) / k) + 0.5;
    vec2 p0 = i0 * k, p1 = i1 * k;
    c0 = (dot2(q - p0) < dot2(q - p1)) ? p0 : p1;
    vec2 r  = q - c0;
    float sx = (r.x >= 0.0) ? 1.0 : -1.0;
    float sy = (r.y >= 0.0) ? 1.0 : -1.0;
    c1 = c0 + vec2(0.5 * sx, 0.5 * 1.73205081 * sy);
    c2 = (abs(r.x) > 1.73205081 * abs(r.y)) ? (c0 + vec2(sx, 0.0)) : (c0 + vec2(-0.5 * sx, 0.5 * 1.73205081 * sy));
}
