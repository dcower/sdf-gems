// ============================================================================
// sdSphericalSegmentSym — Exact 3D Two-Plane Cut Sphere (Conjugate-Square Reduction)
// @author Discovered using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`, 0/1/2 `sqrt`s down from 2/3 `sqrt`s, `0.0` error)
// ============================================================================

float sdSphericalSegmentSym( vec3 p, float r, float h )
{
    float x2 = dot(p.xz, p.xz), y = abs(p.y), w2 = r*r - h*h, d = x2 - w2, dy = y - h;
    if( w2*y*y < h*h*x2 && (r - h)*(r + h + 2.0*dy) < x2 )
        return sqrt(x2 + y*y) - r;
    return (d < 0.0) ? dy : sqrt(d*d / (x2 + w2 + 2.0*sqrt(w2*x2)) + dy*dy);
}
