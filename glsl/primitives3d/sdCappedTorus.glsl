// ============================================================================
// sdCappedTorus — Exact 3D Capped Torus Arc
// @author Inigo Quilez
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes
// ============================================================================

float sdCappedTorus( vec3 p, vec2 sc, float ra, float rb )
{
    p.x = abs(p.x);
    float k = (sc.y*p.x > sc.x*p.y) ? dot(p.xy, sc) : length(p.xy);
    return sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
}
