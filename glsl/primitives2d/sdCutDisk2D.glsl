// ============================================================================
// sdCutDisk2D — 0/1-Sqrt 2D Cut Disk, Cut Hollow Disk, & Analytical Gradient
// @author Inigo Quilez — upgraded via Signed-Square Reduction using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1` everywhere; 0 sqrts in chord strip, 1 sqrt in disk region; down from 1-2 sqrts)
// ============================================================================
//
// Mathematical Foundation:
//   In Inigo Quilez's baseline `sdCutDisk(vec2 p, float r, float h)`
//   (`distfunctions2d.txt` lines 1036-1069), the chord half-width
//     float w = sqrt(r*r - h*h);   // sqrt #1
//   is computed unconditionally on every pixel before branching into:
//     float s = max( (h - r)*p.x*p.x + w*w*(h + r - 2.0*p.y), h*p.x - w*p.y );
//     return (s < 0.0) ? length(p) - r : (p.x < w) ? h - p.y : length(p - vec2(w, h));
//
//   1. Signed-Square & Conjugate-Square Reduction (0 Sqrts in Chord Strip, 1 Sqrt in Disk):
//      Unlike 3D `sdCutSphere` (where radial distance x = length(p.xz) needed a
//      sqrt), in 2D p.x = abs(p.x) is already rational!
//      Let w2 = r*r - h*h and x2 = p.x*p.x. By the Signed-Square Monotonicity
//      Identity:
//        h * p.x < w * p.y   <==>   h * abs(h) * x2 < w2 * p.y * abs(p.y),
//        p.x < w             <==>   x2 < w2,
//      and by Conjugate-Square factorization:
//        (h - r)*x2 + w2*(h + r - 2*p.y) < 0   <==>   (r + h)*(r - h + 2*(h - p.y)) < x2.
//      This completely eliminates `w = sqrt(r*r - h*h)` from the region tests and
//      from the first two return branches:
//        - Flat chord normal strip (x2 < w2):  0 sqrts (returns h - p.y, down from 1 sqrt!)
//        - Circular disk Voronoi region:       1 sqrt  (returns sqrt(x2 + p.y^2) - r, down from 2 sqrts!)
//        - Corner wedge (x2 >= w2):            localized sqrt(w2) only inside the corner wedge!

// 1. 0/1-Sqrt 2D Cut Disk (0 sqrts in chord strip, 1 sqrt in circular disk region)
float sdCutDisk2D( in vec2 p, in float r, in float h )
{
    float w2 = r * r - h * h;
    p.x = abs(p.x);
    float x2 = p.x * p.x, dy = h - p.y;
    if( (r + h) * (r - h + 2.0 * dy) < x2 && h * abs(h) * x2 < w2 * p.y * abs(p.y) )
        return sqrt(x2 + p.y * p.y) - r; // 1 sqrt (down from 2)!
    return (x2 < w2) ? dy : length(vec2(p.x - sqrt(w2), dy)); // 0 sqrts in chord strip!
}

// 2. 1-Sqrt-Shell 2D Cut Hollow Disk (eliminates unconditional w = sqrt(r*r - h*h) in circular shell)
float sdCutHollowDisk2D( in vec2 p, in float r, in float h, in float t )
{
    float w2 = r * r - h * h;
    p.x = abs(p.x);
    float x2 = p.x * p.x;
    return ((h * abs(h) * x2 < w2 * p.y * abs(p.y))
        ? length(vec2(p.x - sqrt(w2), p.y - h))
        : abs(sqrt(x2 + p.y * p.y) - r)) - t;
}

// 3. Analytical 2D Distance and Unit Gradient vec3(d, grad.xy) (0/1-SFU in chord strip & disk)
vec3 sdgCutDisk2D( in vec2 p, in float r, in float h )
{
    float sx = (p.x < 0.0) ? -1.0 : 1.0;
    p.x = abs(p.x);
    float w2 = r * r - h * h;
    float x2 = p.x * p.x, dy = h - p.y;
    if( (r + h) * (r - h + 2.0 * dy) < x2 && h * abs(h) * x2 < w2 * p.y * abs(p.y) )
    {
        float r2   = x2 + p.y * p.y;
        float invR = inversesqrt(max(r2, 1e-30));
        return vec3(r2 * invR - r, sx * p.x * invR, p.y * invR); // 1 inversesqrt!
    }
    if( x2 < w2 )
        return vec3(dy, 0.0, -1.0); // 0 SFU ops in chord strip!
    vec2  v    = vec2(p.x - sqrt(w2), -dy);
    float v2   = dot(v, v);
    float invV = inversesqrt(max(v2, 1e-30));
    return vec3(v2 * invV, sx * v.x * invV, v.y * invV);
}
