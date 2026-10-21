// ============================================================================
// sdCutHollowSphere — Exact 3D Cut Hollow Spherical Shell (Conjugate-Square & Signed-Square Reduction)
// @author Inigo Quilez — upgraded/optimized using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions/
// @exact  Yes (`||grad d|| = 1`, 1/2 `sqrt`s down from 3 `sqrt`s, `0.0` error)
// ============================================================================

float sdCutHollowSphere( vec3 p, float r, float h, float t )
{
    float x2 = dot(p.xz, p.xz), w2 = r*r - h*h, d = x2 - w2, dy = p.y - h;
    return ((h*abs(h)*x2 < w2*p.y*abs(p.y))
        ? sqrt(d*d / max(x2 + w2 + 2.0*sqrt(w2*x2), 1e-30) + dy*dy)
        : abs(sqrt(x2 + p.y*p.y) - r)) - t;
}
