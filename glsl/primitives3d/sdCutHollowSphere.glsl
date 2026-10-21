// ============================================================================
// sdCutHollowSphere — Exact 3D Cut Hollow Spherical Shell
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdCutHollowSphere( vec3 p, float r, float h, float t )
{
    float w = sqrt(r*r - h*h);
    vec2 q = vec2( length(p.xz), p.y );
    return ((h*q.x < w*q.y) ? length(q - vec2(w,h)) : abs(length(q) - r)) - t;
}
