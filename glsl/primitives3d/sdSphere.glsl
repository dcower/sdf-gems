// ============================================================================
// sdSphere — Exact 3D Sphere
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdSphere( vec3 p, float s )
{
    return length(p) - s;
}
