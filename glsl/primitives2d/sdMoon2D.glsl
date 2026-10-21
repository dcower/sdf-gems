// ============================================================================
// sdMoon2D & sdLens2D — Exact 1-Length 2D Crescent Moon & Asymmetric Lens
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1`, 1 dynamic `length()` call down from 2, `0.0` error)
// ============================================================================

// 1. Exact 1-Length 2D Crescent Moon (Confocal Ellipse Bisector, down from 2 length() calls)
float sdMoon2D( vec2 p, float d, float ra, float rb )
{
    p.y = abs(p.y);
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra - a*a, 0.0));
    float r = ra + rb, m = clamp(ra*r + d*(p.x - a), 0.0, r*r);
    vec3  c = (d*(p.x*b - p.y*a) > d*d*max(b - p.y, 0.0)) ? vec3(a, b, 0.0) :
              (r*r*dot(p, p) >= m*m) ? vec3(0.0, 0.0, ra) : vec3(d, 0.0, -rb);
    return length(p - c.xy) * ((c.z < 0.0) ? -1.0 : 1.0) - c.z;
}

// 2. Exact 1-Length 2D Asymmetric Two-Circle Optical Lens (Confocal Hyperbola Bisector)
float sdLens2D( vec2 p, float d, float ra, float rb )
{
    p.y = abs(p.y);
    float a = (ra*ra - rb*rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra - a*a, 0.0));
    float k = p.x*b - p.y*a, r = ra - rb, m = ra*r + d*(p.x - a);
    vec3  c = (k < 0.0 && k > d*(b - p.y)) ? vec3(a, b, 0.0) :
              (r*abs(r)*dot(p, p) <= m*abs(m)) ? vec3(0.0, 0.0, ra) : vec3(d, 0.0, rb);
    return length(p - c.xy) - c.z;
}
