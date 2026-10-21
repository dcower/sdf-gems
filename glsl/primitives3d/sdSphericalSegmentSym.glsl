// ============================================================================
// sdSphericalSegmentSym — Exact 3D Two-Plane Cut Sphere (Spherical Segment)
// @author Dillon Cower (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (`||grad d|| = 1`)
// ============================================================================

float sdSphericalSegmentSym( vec3 p, float r, float h )
{
    float w = sqrt(r*r - h*h);
    vec2  q = vec2( length(p.xz), abs(p.y) );
    float s = max( w*q.y - h*q.x, (r - h)*(r - h + 2.0*q.y) - q.x*q.x );
    return (s < 0.0) ? length(q) - r :
           (q.x < w) ? q.y - h       :
                       length(q - vec2(w, h));
}
