// ============================================================================
// sdDeathStar — Exact 2-Length 3D Sphere Subtracted by Sphere (Confocal Ellipse Bisector)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 2 dynamic `length()` calls down from 3, `0.0` error)
// ============================================================================

float sdDeathStar( vec3 p2, float ra, float rb, float d )
{
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra - a*a, 0.0));
    vec2  p = vec2( p2.x, length(p2.yz) );
    float r = ra + rb, m = clamp(ra*r + d*(p.x - a), 0.0, r*r);
    vec3  c = (p.x*b - p.y*a > d*max(b - p.y, 0.0)) ? vec3(a, b, 0.0) :
              (r*r*dot(p, p) >= m*m) ? vec3(0.0, 0.0, ra) : vec3(d, 0.0, -rb);
    return length(p - c.xy) * ((c.z < 0.0) ? -1.0 : 1.0) - c.z;
}
