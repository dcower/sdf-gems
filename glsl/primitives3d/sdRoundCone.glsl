// ============================================================================
// sdRoundCone — Exact 3D Round Cone (Swept Sphere-Tapered Capsule)
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdRoundCone( vec3 p, float r1, float r2, float h )
{
    float b = (r1 - r2)/h;
    float a = sqrt(1.0 - b*b);
    vec2 q = vec2( length(p.xz), p.y );
    float k = dot(q, vec2(-b, a));
    if( k < 0.0 ) return length(q) - r1;
    if( k > a*h ) return length(q - vec2(0.0, h)) - r2;
    return dot(q, vec2(a, b)) - r1;
}
