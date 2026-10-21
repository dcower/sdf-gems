// ============================================================================
// opElongate — Exact 3D Domain Elongation
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

vec4 opElongate( in vec3 p, in vec3 h )
{
    vec3 q = abs(p) - h;
    return vec4( max(q, 0.0), min(max(q.x, max(q.y, q.z)), 0.0) );
}
