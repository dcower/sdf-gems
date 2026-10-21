// ============================================================================
// sdLens — Exact 2-Length 3D Asymmetric Two-Sphere Optical Lens (Confocal Hyperbola Bisector)
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, 2 dynamic `length()` calls down from 3, `0.0` error)
// ============================================================================

float sdLens( vec3 p2, float ra, float rb, float d )
{
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra - a*a, 0.0));
    vec2  p = vec2( p2.x, length(p2.yz) );
    float k = p.x*b - p.y*a, r = ra - rb, m = ra*r + d*(p.x - a);
    vec3  c = (k < 0.0 && k > d*(b - p.y)) ? vec3(a, b, 0.0) :
              (r*abs(r)*dot(p, p) <= m*abs(m)) ? vec3(0.0, 0.0, ra) : vec3(d, 0.0, rb);
    return length(p - c.xy) - c.z;
}
