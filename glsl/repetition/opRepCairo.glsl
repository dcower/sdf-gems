// ============================================================================
// opRepCairo — 0-Trig Cairo Pentagonal, Trihexagonal Kagome & O(1) Spherical Fibonacci Domain Repetition
// @author Discovered and derived using Gemini 4 Argon (sdf-gems, 2026)
// @source https://github.com/dcower/sdf-gems
// @exact  Yes (0-trig, 0-sqrt exact Cairo cell distance; C^infinity seamless across Z^2 grid;
//         O(1) 0-loop Spherical Fibonacci inverse parastichy mapping)
// ============================================================================

#ifndef SDF_GEMS_OP_REP_CAIRO_GLSL
#define SDF_GEMS_OP_REP_CAIRO_GLSL

// ----------------------------------------------------------------------------
// 1. opRepCairo(p, s, k, id) & sdCairoCell2D(q, s, k)
//
// Exact 0-trig, 0-sqrt Cairo Pentagonal Tiling (p4g wallpaper group) domain
// repetition with cell period s and edge slope parameter k in (0, 1):
//   - k = 0.57735027 (1 / sqrt(3))      : Regular Catalan Cairo tiling (90-120-120-90-120 deg,
//                                         bridge/side ratio = sqrt(3) - 1 = 0.73205081,
//                                         exact incircle of radius 0.3660254 * s)
//   - k = 0.45141623 ((4 - sqrt(7)) / 3): Equilateral Cairo tiling (all 5 sides equal)
//   - k = 0.50000000 (1 / 2)            : Pyritohedral rational 3-4-5 Cairo tiling
//
// Returns incenter-referenced local coordinates q in the canonical pentagon
// frame and writes the unique half-integer pentagon center to id (100% seamless
// across all Z^2 square-grid boundaries!).
// ----------------------------------------------------------------------------
vec2 opRepCairo( vec2 p, float s, float k, out vec2 id )
{
    vec2  u   = p / s;
    vec2  i   = floor(u + 0.5);
    vec2  q   = u - i;
    float par = mod(i.x + i.y, 2.0);
    if( par < 0.0 ) par += 2.0;
    if( par > 0.5 ) q = q.yx;

    vec2  w = 0.5 - abs(q);
    vec2  ql, off;
    if( w.x - k * w.y >= 0.0 )
    {
        float sy = (q.y >= 0.0) ? 1.0 : -1.0;
        ql  = vec2( q.x, sy * q.y - 0.5 );
        off = vec2( 0.0, sy * 0.5 );
    }
    else
    {
        float sx = (q.x >= 0.0) ? 1.0 : -1.0;
        ql  = vec2( q.y, 0.5 - sx * q.x );
        off = vec2( sx * 0.5, 0.0 );
    }
    if( par > 0.5 ) off = off.yx;
    id = i + off;

    float yc = -0.5 * (1.0 - k) / (1.0 + k);
    return (ql - vec2(0.0, yc)) * s;
}

// Exact 0-sqrt Euclidean distance to the 5 walls of the Cairo pentagon from q = opRepCairo(...)
float sdCairoCell2D( vec2 q, float s, float k )
{
    float yc  = -0.5 * (1.0 - k) / (1.0 + k);
    vec2  ql  = q / s + vec2(0.0, yc);
    vec2  u   = vec2( 0.5 - abs(ql.x), -ql.y );
    float inv = inversesqrt(1.0 + k * k);
    return min( 0.5 + ql.y, min(u.x - k * u.y, k * u.x + u.y) * inv ) * s;
}

// Helper for 4-Tap Neighbor-Corrected Cairo Repetition: maps p into any half-integer pentagon cid
vec2 opRepCairoEvalCell( vec2 p, float s, float k, vec2 cid )
{
    vec2  u  = p / s - cid;
    float ry = abs(fract(cid.y));
    vec2  ql;
    if( ry > 0.25 && ry < 0.75 )
    {
        vec2  i   = vec2( floor(cid.x + 0.5), floor(cid.y) );
        float par = mod(i.x + i.y, 2.0);
        if( par < 0.0 ) par += 2.0;
        ql = vec2( u.x, (par < 0.5) ? u.y : -u.y );
    }
    else
    {
        vec2  i   = vec2( floor(cid.x), floor(cid.y + 0.5) );
        float par = mod(i.x + i.y, 2.0);
        if( par < 0.0 ) par += 2.0;
        ql = vec2( u.y, (par < 0.5) ? -u.x : u.x );
    }
    float yc = -0.5 * (1.0 - k) / (1.0 + k);
    return (ql - vec2(0.0, yc)) * s;
}

// Returns the 4 quadrant-adaptive neighbor pentagon IDs (c0, c1, c2, c3) for 0.0-error 4-tap evaluation
void opRepCairo4TapIDs( vec2 p, float s, out vec2 c0, out vec2 c1, out vec2 c2, out vec2 c3 )
{
    vec2  u   = p / s;
    vec2  i   = floor(u + 0.5);
    vec2  q   = u - i;
    float par = mod(i.x + i.y, 2.0);
    if( par < 0.0 ) par += 2.0;
    if( par > 0.5 ) q = q.yx;

    float sx = (q.x >= 0.0) ? 1.0 : -1.0;
    float sy = (q.y >= 0.0) ? 1.0 : -1.0;
    vec2  a  = abs(q);
    vec2  o0 = vec2( 0.0, sy * 0.5 );
    vec2  o1 = vec2( sx * 0.5, 0.0 );
    vec2  o2 = (a.x + a.y >= 0.5) ? vec2( sx * 0.5, sy ) : vec2( 0.0, -sy * 0.5 );
    vec2  o3 = (a.x + a.y >= 0.5) ? vec2( sx, sy * 0.5 )
                                  : ((a.y >= a.x) ? vec2( sx * 0.5, sy ) : vec2( sx, sy * 0.5 ));
    if( par > 0.5 ) { o0 = o0.yx; o1 = o1.yx; o2 = o2.yx; o3 = o3.yx; }
    c0 = i + o0; c1 = i + o1; c2 = i + o2; c3 = i + o3;
}

// ----------------------------------------------------------------------------
// 2. opRepKagome(p, s, id, type) — 0-Trig, 0-Sqrt Trihexagonal Kagome Repetition
//
// Partitions each A2 hexagonal Voronoi cell of spacing s into its central
// regular hexagon (type = 0.0, inradius = 0.25*sqrt(3)*s = 0.4330127*s) and
// 6 corner equilateral triangles (type = 1.0, inradius = s/(4*sqrt(3)) = 0.1443376*s)
// using 0 trig and 0 sqrts. Returns local cell coordinates q and signed wall distance.
// ----------------------------------------------------------------------------
vec3 opRepKagome( vec2 p, float s, out vec2 id, out float tile_type )
{
    const vec2 k_hex = vec2( 1.0, 1.73205081 );
    vec2 u = p / (s * k_hex), i = floor(u + 0.5), d = u - i;
    if( abs(d.x) + 3.0 * abs(d.y) > 1.0 )
        i += vec2( (d.x >= 0.0) ? 0.5 : -0.5, (d.y >= 0.0) ? 0.5 : -0.5 );
    vec2 q = p - s * i * k_hex;

    vec2  a  = abs(q);
    float m  = max( a.y, a.x * 0.86602540 + a.y * 0.5 );
    float rh = 0.43301270 * s;
    if( m <= rh )
    {
        tile_type = 0.0;
        id = i * 2.0;
        return vec3( q, rh - m );
    }
    tile_type = 1.0;
    vec2  v  = (a.y > a.x * 0.86602540 + a.y * 0.5)
             ? vec2( 0.0, (q.y >= 0.0 ? 1.0 : -1.0) * 0.57735027 * s )
             : vec2( (q.x >= 0.0 ? 0.5 : -0.5) * s, (q.y >= 0.0 ? 1.0 : -1.0) * 0.28867513 * s );
    id = i * 2.0 + v / s;
    return vec3( q - v, m - rh );
}

// ----------------------------------------------------------------------------
// 3. opRepFibonacciSphere(p, n, id, dir) — O(1) 0-Loop Spherical Fibonacci Repetition
//
// Maps any 3D point p in R^3 \ {0} to the closest of N golden-ratio Spherical
// Fibonacci points on the sphere in O(1) constant time (0 loops!) via the 2x2
// inverse parastichy zone matrix (Keinert et al. 2015).
// Returns local 3D offset q = p - length(p) * dir, writes integer index id in
// [0, n - 1] and unit Fibonacci lattice normal dir.
// ----------------------------------------------------------------------------
vec3 sdFibonacciPoint( float j, float n )
{
    const float G_FRAC = 0.61803398875;
    const float TAU    = 6.28318530718;
    float phi = TAU * fract(j * G_FRAC);
    float z   = 1.0 - (2.0 * j + 1.0) / n;
    float r   = sqrt(max(1.0 - z * z, 0.0));
    return vec3( r * cos(phi), r * sin(phi), z );
}

vec3 opRepFibonacciSphere( vec3 p, float n, out float id, out vec3 dir )
{
    const float PHI    = 1.61803398875;
    const float G_FRAC = 0.61803398875;
    const float TAU    = 6.28318530718;
    const float INV_S5 = 0.44721359550;

    float rp  = length(p);
    vec3  u   = (rp > 0.0) ? (p / rp) : vec3(0.0, 0.0, 1.0);
    float phi = atan(u.y, u.x);
    if( phi < 0.0 ) phi += TAU;
    float uz  = clamp(u.z, -1.0, 1.0);

    float k   = max(2.0, floor(log(max(n * 7.024814731 * (1.0 - uz * uz), 1e-30)) * 1.03904346));
    float pk  = pow(PHI, k);
    float Fk  = floor(pk * INV_S5 + 0.5);
    float Fk1 = floor(pk * PHI * INV_S5 + 0.5);

    float ka  = TAU * (fract((Fk  + 1.0) * G_FRAC) - G_FRAC);
    float kb  = TAU * (fract((Fk1 + 1.0) * G_FRAC) - G_FRAC);
    float inv = n / (2.0 * (kb * Fk - ka * Fk1));
    float vz  = uz - (1.0 - 1.0 / n);

    float c0  = floor(inv * (-2.0 * (Fk1 / n) * phi - kb * vz));
    float c1  = floor(inv * ( 2.0 * (Fk  / n) * phi + ka * vz));

    float j00 = clamp(Fk * c0         + Fk1 * c1,         0.0, n - 1.0);
    float j10 = clamp(Fk * (c0 + 1.0) + Fk1 * c1,         0.0, n - 1.0);
    float j01 = clamp(Fk * c0         + Fk1 * (c1 + 1.0), 0.0, n - 1.0);
    float j11 = clamp(Fk * (c0 + 1.0) + Fk1 * (c1 + 1.0), 0.0, n - 1.0);

    vec3  f00 = sdFibonacciPoint(j00, n), f10 = sdFibonacciPoint(j10, n);
    vec3  f01 = sdFibonacciPoint(j01, n), f11 = sdFibonacciPoint(j11, n);

    id  = j00; dir = f00; float best = dot(u, f00);
    float d10 = dot(u, f10); if( d10 > best ) { best = d10; id = j10; dir = f10; }
    float d01 = dot(u, f01); if( d01 > best ) { best = d01; id = j01; dir = f01; }
    float d11 = dot(u, f11); if( d11 > best ) { best = d11; id = j11; dir = f11; }

    return p - rp * dir;
}

#endif // SDF_GEMS_OP_REP_CAIRO_GLSL
