// ============================================================================
// sdStar2D — 1-Sqrt Orthonormal-Weyl 2D Star, Regular n-Gon, & 3D Regular Prism
// @author Inigo Quilez — upgraded via Orthonormal Weyl Decomposition using Gemini 4 Argon (sdf-gems, 2026)
// @source https://iquilezles.org/articles/distfunctions2d/
// @exact  Yes (`||grad d|| = 1` everywhere; eliminates the 2nd post-polar-fold sqrt across all 2n star edge strips & 100% of n-gon interior/strips)
// ============================================================================
//
// Mathematical Foundation:
//   In Inigo Quilez's baseline `sdStar(vec2 p, float r, int n, float m)`
//   (`distfunctions2d.txt` lines 926-997), `length(p)` (sqrt #1) is used for
//   polar sector folding into `[-pi/n, pi/n]`, followed by clamping onto the
//   edge direction `ecs = vec2(cos(en), sin(en))` and calling `length(p)*sign(p.x)`
//   (sqrt #2) unconditionally at every pixel.
//
//   1. Orthonormal Weyl Decomposition (`sdStar2D`):
//      Since `ecs` is a unit vector, `(-ecs, (ecs.y, -ecs.x))` forms an exact
//      orthonormal basis! Let `h = 0.5 * r * acs.y / ecs.y` be the half-edge
//      length, with tangential overshoot `d.x = abs(-dot(p, ecs) - h) - h` and
//      signed outward normal distance `d.y = p.x * ecs.y - p.y * ecs.x`.
//      Whenever `d.x <= 0.0` (across all 2*n edge-normal strips of the star),
//      the closest point lies strictly in the edge interior, so the exact signed
//      distance is simply `d.y` — eliminating the second `sqrt` completely!
//
//   2. Regular n-Gon (`sdRegularPolygon2D`) & 1-Sqrt 3D `sdRegularPrism`:
//      For a regular n-gon (`ecs = vec2(0, 1)`), the folded edge is vertical at
//      `p.x = 0` with `p.y in [-r*sin(pi/n), 0]`. Because the chamber wedge is
//      acute (`pi/n < pi/2`), `p.y > 0.0` implies `p.x > 0.0` strictly, so
//      100% of the polygon interior and all n exterior edge strips have
//      `p.y <= 0.0` and return `p.x` with ZERO post-fold sqrts (and at most 1
//      post-fold sqrt in 3D `sdRegularPrism`)!

// 1. 1-Sqrt-Strip 2D Star Polygon (n points, angle divisor m in [2, n])
float sdStar2D( in vec2 p, in float r, in int n, in float m )
{
    float an  = 3.14159265359 / float(n);
    float en  = 3.14159265359 / m;
    vec2  acs = vec2(cos(an), sin(an));
    vec2  ecs = vec2(cos(en), sin(en));
    float bn  = mod(atan(p.x, p.y), 2.0 * an) - an;
    p = length(p) * vec2(cos(bn), abs(sin(bn))) - r * acs;
    float h = 0.5 * r * acs.y / ecs.y;
    vec2  d = vec2(abs(-dot(p, ecs) - h) - h, p.x * ecs.y - p.y * ecs.x);
    return (d.x > 0.0) ? length(d) * ((d.y < 0.0) ? -1.0 : 1.0) : d.y; // 0 post-fold sqrts in all 2n edge strips!
}

// 2. 1-Sqrt-Interior 2D Regular n-Gon (circumradius r, n >= 3 sides)
float sdRegularPolygon2D( in vec2 p, in float r, in int n )
{
    float an  = 3.14159265359 / float(n);
    vec2  acs = vec2(cos(an), sin(an));
    float bn  = mod(atan(p.x, p.y), 2.0 * an) - an;
    p = length(p) * vec2(cos(bn), abs(sin(bn))) - r * acs;
    return (p.y > 0.0) ? length(p) : p.x; // 0 post-fold sqrts in 100% of interior & n edge strips!
}

// 3. Exact 3D Regular n-Gon Prism (circumradius r, half-height hz, n >= 3 sides)
float sdRegularPrism( in vec3 p, in float r, in float hz, in int n )
{
    float an  = 3.14159265359 / float(n);
    vec2  acs = vec2(cos(an), sin(an));
    float bn  = mod(atan(p.x, p.y), 2.0 * an) - an;
    vec2  q   = length(p.xy) * vec2(cos(bn), abs(sin(bn))) - r * acs;
    vec3  d   = vec3(q.x, max(q.y, 0.0), abs(p.z) - hz);
    return min(max(d.x, d.z), 0.0) + length(max(d, 0.0)); // At most 1 post-fold sqrt in R^3!
}
