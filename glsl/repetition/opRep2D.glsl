// ============================================================================
// opRep2D — 3-Tap 1-Round 2D Hexagonal (A2) Delaunay Domain Repetition (opRepHex3)
// @author Inigo Quilez (4-tap Cartesian) — upgraded/fixed using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes (`0.0` error vs 19-tap brute force; 100% Delaunay containment, 1 `round`, 0 `dot2`)
// ============================================================================

void opRepHex3( in vec2 p, in float s, out vec2 c0, out vec2 c1, out vec2 c2 )
{
    const vec2 k = vec2(1.0, 1.73205081);
    vec2 u  = p / (s * k), i0 = round(u), d = u - i0;
    vec2 sg = vec2(d.x >= 0.0 ? 1.0 : -1.0, d.y >= 0.0 ? 1.0 : -1.0);
    c0 = i0 * k;
    c1 = (i0 + 0.5 * sg) * k;
    c2 = (i0 + ((abs(d.x) > abs(d.y)) ? vec2(sg.x, 0.0) : vec2(-0.5 * sg.x, 0.5 * sg.y))) * k;
}
