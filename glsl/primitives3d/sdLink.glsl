// ============================================================================
// sdLink — Exact 3D Chain Link
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdLink( vec3 p, float le, float r1, float r2 )
{
    vec3 q = vec3( p.x, max(abs(p.y) - le, 0.0), p.z );
    return length(vec2(length(q.xy) - r1, q.z)) - r2;
}
