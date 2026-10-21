// ============================================================================
// sdSphericalZone — Exact 3D Hollow Spherical Zone (Conjugate-Square & Signed-Square Reduction)
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, 1/2 `sqrt`s down from 4 `sqrt`s, `0.0` error)
// ============================================================================

float sdSphericalZone( vec3 p, float r, float h1, float h2, float t )
{
    float x2 = dot(p.xz, p.xz), ys = p.y*abs(p.y);
    float w1 = r*r - h1*h1, w2 = r*r - h2*h2;
    bool  c1 = h1*abs(h1)*x2 > w1*ys;
    bool  c2 = !c1 && (w2*ys > h2*abs(h2)*x2);
    if( !c1 && !c2 ) return abs(sqrt(x2 + p.y*p.y) - r) - t;
    float w = c1 ? w1 : w2, h = c1 ? h1 : h2, d = x2 - w, dy = p.y - h;
    return sqrt(d*d / max(x2 + w + 2.0*sqrt(w*x2), 1e-30) + dy*dy) - t;
}
