// ============================================================================
// sdCutSphere — Exact 3D Single-Plane Cut Sphere
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdCutSphere( vec3 p, float r, float h )
{
    float w = sqrt(r*r - h*h);
    vec2 q = vec2( length(p.xz), p.y );
    float s = max( (h - r)*q.x*q.x + w*w*(h + r - 2.0*q.y), h*q.x - w*q.y );
    return (s < 0.0) ? length(q) - r :
           (q.x < w) ? h - q.y       :
                       length(q - vec2(w,h));
}
