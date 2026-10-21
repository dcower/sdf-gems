// ============================================================================
// sdCone — Exact 3D Cone
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdCone( vec3 p, vec2 c, float h )
{
    vec2 q = h*vec2(c.x/c.y, -1.0);
    vec2 w = vec2( length(p.xz), p.y );
    vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
    vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
    float k = sign( q.y );
    float d = min(dot( a, a ), dot(b, b));
    float s = max( k*(w.x*q.y - w.y*q.x), k*(w.y - q.y) );
    return sqrt(d)*sign(s);
}
