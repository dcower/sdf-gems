// ============================================================================
// sdCircle — Exact 2D Circle
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdCircle( vec2 p, float r )
{
    return length(p) - r;
}
