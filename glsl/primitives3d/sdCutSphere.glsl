// ============================================================================
// sdCutSphere — Exact 3D Single-Plane Cut Sphere (Conjugate-Square & Signed-Square Reduction)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 0/1/2 `sqrt`s down from 2/3 `sqrt`s, `0.0` error)
// ============================================================================

float sdCutSphere( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), w2 = r*r - h*h, d = x2 - w2, dy = h - p.y;
    if( (r + h)*(r - h + 2.0*dy) < x2 && h*abs(h)*x2 < w2*p.y*abs(p.y) )
        return sqrt(x2 + p.y*p.y) - r;
    return (d < 0.0) ? dy : sqrt(d*d / (x2 + w2 + 2.0*sqrt(w2*x2)) + dy*dy);
}
