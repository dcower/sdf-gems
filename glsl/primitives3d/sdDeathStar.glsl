// ============================================================================
// sdDeathStar — Exact 3D Sphere Subtracted by Sphere (Death Star)
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdDeathStar( vec3 p2, float ra, float rb, float d )
{
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra - a*a, 0.0));
    vec2 p = vec2( p2.x, length(p2.yz) );
    if( p.x*b - p.y*a > d*max(b - p.y, 0.0) )
        return length(p - vec2(a,b));
    else
        return max( (length(p) - ra), -(length(p - vec2(d,0.0)) - rb));
}
