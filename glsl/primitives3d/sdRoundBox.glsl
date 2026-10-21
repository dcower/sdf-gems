// ============================================================================
// sdRoundBox — Exact 3D Rounded Box
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdRoundBox( vec3 p, vec3 b, float r )
{
    vec3 q = abs(p) - b + r;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0) - r;
}
