// ============================================================================
// opRep2D — 2D Neighbor-Corrected Domain Repetition (4-Tap Cartesian)
// @author Inigo Quilez
// @source https://iquilezles.org/articles/sdfrepetition/
// @exact  Yes
// ============================================================================

void opRep2D_4Tap( in vec2 p, in float s, out vec2 c0, out vec2 c1, out vec2 c2, out vec2 c3 )
{
    vec2 q = p / s;
    vec2 id = round(q);
    vec2 o = sign(q - id);
    c0 = id + vec2(0.0, 0.0);
    c1 = id + vec2(o.x, 0.0);
    c2 = id + vec2(0.0, o.y);
    c3 = id + vec2(o.x, o.y);
}
