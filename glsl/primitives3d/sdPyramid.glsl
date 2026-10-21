// ============================================================================
// sdPyramid — Exact 3D Square Pyramid (Base-Face Bug Fixed)
// @author Inigo Quilez — upgraded/fixed by Dillon Cower (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1` everywhere in R^3)
// ============================================================================

float sdPyramid( vec3 p, float h )
{
    float m2 = h*h + 0.25;
    p.xz = abs(p.xz);
    p.xz = (p.z > p.x) ? p.zx : p.xz;
    p.xz -= 0.5;

    vec3 q = vec3( p.z, h*p.y - 0.5*p.x, h*p.x + 0.5*p.y );
    float s = max(-q.x, 0.0);
    float t = clamp( (q.y - 0.5*p.z)/(m2 + 0.25), 0.0, 1.0 );

    float a = m2*(q.x + s)*(q.x + s) + q.y*q.y;
    float b = m2*(q.x + 0.5*t)*(q.x + 0.5*t) + (q.y - m2*t)*(q.y - m2*t);
    float d2 = min(q.y, -q.x*m2 - q.y*0.5) > 0.0 ? 0.0 : min(a, b);

    float d = (d2 + q.z*q.z) / m2;
    return sqrt( max(q.z, p.x) < 0.0 ? min(d, p.y*p.y) : d ) * sign(max(q.z, -p.y));
}
