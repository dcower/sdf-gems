// ============================================================================
// ptHarmonic — Closed-Form Capacitary (Equilibrium) Potentials of Primitives
//   ptSphere, ptSpheroid, ptDisk, ptSegment, ptEllipsoid (+ ptEllipsoidCapacity,
//   ptEllipsoidCap, ptCarlsonRF), ptTwoSpheres, ptConvexBound, ptCapacityFit,
//   ptUnionCoupled
//
// u(p) = probability that Brownian motion started at p hits the body K before
// escaping to infinity = the harmonic function outside K with u = 1 on K and
// u -> 0 at infinity (u ~ C / |p|, C = capacity of K in units of length).
// Functions returning vec4 give (u, grad u). Inside K they return (1, 0, 0, 0).
//
// @author Claude Fable 5.2 (sdf-gems, 2026). The formulas are classical
//         potential theory (conducting ellipsoid, image charges); see the doc.
// @exact  ptSphere, ptSpheroid, ptDisk, ptSegment (a prolate spheroid, NOT a capsule): Exact
//         potentials (GPU FP32 vs FP64 references: max |du| <= 2.8e-6 on >= 5,700 exterior points
//         each; consistent with GPU Walk-on-Spheres on the repo SDFs).
//         ptEllipsoid: Exact formula with a fixed 6-step Newton solve. Verified (max |du| <= 7.2e-7)
//         for axis ratios up to 10:1. Approximate beyond: 6.8e-5 at 33:1, 5.7e-4 at 50:1,
//         1.5e-2 at 100:1 (errors concentrated next to the thin edge).
//         ptTwoSpheres: Approximate (image series truncated at PT_TWO_SPHERES_IMAGES = 12 terms):
//         max |du| 1.1e-5 at gap = 1 radius, 7.3e-3 at gap = 0.1, 2.5e-2 touching; the result can
//         exceed 1 by the same amount near the gap (not clamped).
//         ptConvexBound: Bounds. No violation in 24,000 samples on 4 convex bodies; the inscribed /
//         circumscribed sphere terms follow from domain monotonicity, the d-dependent upper term
//         is only tested, not proved here.
//         ptCapacityFit, ptUnionCoupled: Approximate (max errors 0.07 - 0.22 and 0.08; see the doc).
//         None of these is a distance: do not sphere-trace them without a gradient bound.
// ============================================================================

#ifndef SDF_GEMS_HAS_PT_HARMONIC
#define SDF_GEMS_HAS_PT_HARMONIC

#ifndef PT_TWO_SPHERES_IMAGES
#define PT_TWO_SPHERES_IMAGES 12
#endif

// atanh(x) / x for 0 <= x < 1 and atan(x) / x for x >= 0, accurate as x -> 0.
float ptAtanhc( float x )
{
    float x2 = x * x;
    return ( x < 0.15 ) ? 1.0 + x2 * (1.0 / 3.0 + x2 * (0.2 + x2 * (1.0 / 7.0 + x2 * (1.0 / 9.0))))
                        : 0.5 * log( (1.0 + x) / (1.0 - x) ) / x;
}
float ptAtanc( float x )
{
    float x2 = x * x;
    return ( x < 0.15 ) ? 1.0 - x2 * (1.0 / 3.0 - x2 * (0.2 - x2 * (1.0 / 7.0 - x2 * (1.0 / 9.0))))
                        : atan( x ) / x;
}

// ---------------------------------------------------------------------------
// Sphere of radius R at the origin: u = R / r. Capacity R.
// ---------------------------------------------------------------------------
vec4 ptSphere( vec3 p, float R )
{
    float r2 = dot( p, p );
    if( r2 <= R * R ) return vec4( 1.0, 0.0, 0.0, 0.0 );
    float ir = inversesqrt( r2 );
    float u  = R * ir;
    return vec4( u, -(u * ir * ir) * p );
}

// ---------------------------------------------------------------------------
// Spheroid of revolution about the y axis: equatorial radius a (xz), polar
// radius c (y); same convention as sdSpheroid(p, vec2(a, c)). c = 0 is the disk.
// lambda = confocal coordinate (largest root of a quadratic, cancellation-free form).
// Capacity: f / atanh(f / c) (prolate), f / atan(f / c) (oblate), f^2 = |a^2 - c^2|.
// ---------------------------------------------------------------------------
vec4 ptSpheroid( vec3 p, float a, float c )
{
    float rho2 = dot( p.xz, p.xz ), y2 = p.y * p.y;
    float a2 = a * a, c2 = c * c;
    float s  = (a - c) * (a + c);                       // > 0 oblate, < 0 prolate
    float b  = rho2 + y2 - a2 - c2;
    float q  = c2 * rho2 + a2 * y2 - a2 * c2;           // > 0 outside
    float w  = rho2 - y2 - s;
    float sd = sqrt( w * w + 4.0 * rho2 * y2 );         // sqrt(b^2 + 4 q), a sum of squares
    float lam = ( b >= 0.0 ) ? 0.5 * (b + sd) : 2.0 * q / (sd - b);
    if( !(lam > 0.0) ) return vec4( 1.0, 0.0, 0.0, 0.0 );

    float A = a2 + lam, C = c2 + lam, sc = sqrt( C );
    float f = sqrt( abs( s ) );
    float u, cap;
    if( s > 0.0 )
    {
        if( f < 0.15 * c ) { cap = c / ptAtanc( f / c ); u = cap * ptAtanc( f / sc ) / sc; }
        else               { cap = f / atan( f, c );     u = atan( f, sc ) / atan( f, c ); }
    }
    else
    {
        cap = c / ptAtanhc( f / c );
        u   = cap * ptAtanhc( f / sc ) / sc;
    }
    // grad u = -(cap / 2) grad(lambda) / ((a^2 + lambda) sqrt(c^2 + lambda))
    vec3  g = vec3( p.x / A, p.y / C, p.z / A );
    float T = rho2 / (A * A) + y2 / (C * C);
    return vec4( u, -(cap / (A * sc * T)) * g );
}

// Thin disk of radius R in the xz plane (normal y): u = (2 / pi) atan(R / sqrt(lambda)). Capacity 2 R / pi.
vec4 ptDisk( vec3 p, float R )
{
    return ptSpheroid( p, R, 0.0 );
}

// ---------------------------------------------------------------------------
// "Needle": prolate spheroid with foci pa, pb and equatorial radius r > 0
// (semi-major axis sqrt(|pb - pa|^2 / 4 + r^2)). The zero-radius segment itself has
// zero capacity in 3D (u = 0 off the segment), so r > 0 is required.
// u = atanh(f / m) / atanh(f / c), m = (|p - pa| + |p - pb|) / 2.
// ---------------------------------------------------------------------------
vec4 ptSegment( vec3 p, vec3 pa, vec3 pb, float r )
{
    vec3  da = p - pa, db = p - pb;
    float la = length( da ), lb = length( db );
    float f  = 0.5 * length( pb - pa );
    float c  = sqrt( f * f + r * r );
    float m  = 0.5 * (la + lb);
    if( m <= c ) return vec4( 1.0, 0.0, 0.0, 0.0 );
    float cap = c / ptAtanhc( f / c );
    float u   = cap * ptAtanhc( f / m ) / m;
    vec3  gm  = 0.5 * (da / la + db / lb);              // la, lb >= c - f > 0 outside
    return vec4( u, -(cap / ((m - f) * (m + f))) * gm );
}

// ---------------------------------------------------------------------------
// Carlson symmetric elliptic integral R_F(x, y, z) (duplication; x, y, z >= 0, at most one zero).
// ---------------------------------------------------------------------------
float ptCarlsonRF( float x, float y, float z )
{
    for( int i = 0; i < 5; i++ )
    {
        float sx = sqrt( x ), sy = sqrt( y ), sz = sqrt( z );
        float l  = sx * (sy + sz) + sy * sz;
        x = 0.25 * (x + l); y = 0.25 * (y + l); z = 0.25 * (z + l);
    }
    float m  = (x + y + z) / 3.0;
    float X  = 1.0 - x / m, Y = 1.0 - y / m, Z = -(X + Y);
    float E2 = X * Y - Z * Z, E3 = X * Y * Z;
    return (1.0 - E2 / 10.0 + E3 / 14.0 + E2 * E2 / 24.0 - 3.0 * E2 * E3 / 44.0) * inversesqrt( m );
}

// Capacity of the ellipsoid with semi-axes r: C = 1 / R_F(rx^2, ry^2, rz^2). Constant per shape: hoist it.
float ptEllipsoidCapacity( vec3 r )
{
    return 1.0 / ptCarlsonRF( r.x * r.x, r.y * r.y, r.z * r.z );
}

// ---------------------------------------------------------------------------
// Triaxial ellipsoid, semi-axes r, capacity cap = ptEllipsoidCapacity(r):
//   u = cap * R_F(rx^2 + lambda, ry^2 + lambda, rz^2 + lambda),
// lambda = largest root of sum p_i^2 / (r_i^2 + lambda) = 1, found by monotone Newton on
// 1 / S(lambda) - 1 (concave, increasing; one step is exact for a sphere) from a lower bound.
// ---------------------------------------------------------------------------
vec4 ptEllipsoidCap( vec3 p, vec3 r, float cap )
{
    vec3 r2 = r * r, p2 = p * p;
    if( dot( p2, 1.0 / r2 ) <= 1.0 ) return vec4( 1.0, 0.0, 0.0, 0.0 );
    float lam = max( dot( p, p ) - max( r2.x, max( r2.y, r2.z ) ), 0.0 );
    vec3  iv;
    float T;
    for( int i = 0; i < 6; i++ )
    {
        iv = 1.0 / (r2 + lam);
        float S = dot( p2, iv );
        T = dot( p2, iv * iv );
        lam += S * (S - 1.0) / T;
    }
    vec3 e = r2 + lam;
    iv = 1.0 / e;
    T  = dot( p2, iv * iv );
    float u = cap * ptCarlsonRF( e.x, e.y, e.z );
    return vec4( min( u, 1.0 ), -(cap * inversesqrt( e.x * e.y * e.z ) / T) * (p * iv) );
}
vec4 ptEllipsoid( vec3 p, vec3 r )
{
    return ptEllipsoidCap( p, r, ptEllipsoidCapacity( r ) );
}

// ---------------------------------------------------------------------------
// Two disjoint (or touching) spheres at the same potential: centres ca, cb, radii ra, rb,
// |cb - ca| >= ra + rb. Two chains of image charges (Kelvin images), PT_TWO_SPHERES_IMAGES
// charges per chain, the last one at half weight (Euler averaging of the alternating tail).
// Also models one sphere above a REFLECTING ground plane (cb = mirror image of ca), with the same
// truncation error (about 2e-2 for a sphere resting on the plane).
// ---------------------------------------------------------------------------
vec4 ptTwoSpheres( vec3 p, vec3 ca, float ra, vec3 cb, float rb )
{
    vec3  ab = cb - ca;
    float D  = length( ab );
    vec3  n  = ab / D;
    vec3  pa = p - ca, pb = p - cb;
    if( dot( pa, pa ) <= ra * ra || dot( pb, pb ) <= rb * rb ) return vec4( 1.0, 0.0, 0.0, 0.0 );

    // Chain 1 starts with charge ra at ca; chain 2 with charge rb at cb. s = distance of the
    // current charge from the centre of the sphere it sits in, measured towards the other sphere.
    float q1 = ra, s1 = 0.0, q2 = rb, s2 = 0.0;
    vec4  acc = vec4( 0.0 );
    for( int i = 0; i < PT_TWO_SPHERES_IMAGES; i++ )
    {
        bool  inA = ( (i & 1) == 0 );                    // chain 1 charge sits in A on even steps
        float wgt = ( i == PT_TWO_SPHERES_IMAGES - 1 ) ? 0.5 : 1.0;
        vec3  x1 = inA ? ca + s1 * n : cb - s1 * n;
        vec3  x2 = inA ? cb - s2 * n : ca + s2 * n;
        vec3  d1 = p - x1, d2 = p - x2;
        float i1 = inversesqrt( dot( d1, d1 ) ), i2 = inversesqrt( dot( d2, d2 ) );
        float u1 = wgt * q1 * i1, u2 = wgt * q2 * i2;
        acc += vec4( u1 + u2, -(u1 * i1 * i1) * d1 - (u2 * i2 * i2) * d2 );
        // Kelvin image in the other sphere (radius R): q' = -q R / t, at R^2 / t from its centre.
        float R1 = inA ? rb : ra, R2 = inA ? ra : rb;
        float t1 = D - s1, t2 = D - s2;
        q1 = -q1 * R1 / t1; s1 = R1 * R1 / t1;
        q2 = -q2 * R2 / t2; s2 = R2 * R2 / t2;
    }
    return acc;
}

// ---------------------------------------------------------------------------
// Bounds (see the header for what is proved and what is only tested) for a convex body K with exact exterior distance d = sdf(p) >= 0, containing the
// ball B(c, rin) and contained in the ball B(c, rc); l = |p - c|. Returns (lower, upper):
//   lower = rin / l                                 (inscribed sphere, domain monotonicity)
//   upper = min( rc / l, 2 rc / sqrt(4 rc^2 + d^2) ) (circumsphere; supporting half-space + diameter)
// ---------------------------------------------------------------------------
vec2 ptConvexBound( float d, float l, float rin, float rc )
{
    if( d <= 0.0 ) return vec2( 1.0 );
    float D = 2.0 * rc;
    return vec2( min( rin / l, 1.0 ), min( min( rc / l, D * inversesqrt( D * D + d * d ) ), 1.0 ) );
}

// One-tap capacity fit u ~ C / (C + d): exact for the sphere, correct far field C / r for every
// body, Approximate elsewhere (see the doc for measured error maps).
float ptCapacityFit( float d, float C )
{
    return C / (C + max( d, 0.0 ));
}

// ---------------------------------------------------------------------------
// Union of two bodies A, B (Approximate): monopole-coupled superposition
//   u ~ alpha uA + beta uB,  alpha + beta * uBatA = 1,  alpha * uAatB + beta = 1,
// where uBatA = u_B at the centre of A and uAatB = u_A at the centre of B (the mean of a harmonic
// function over a sphere is its value at the centre, so the mean boundary value is exact for spheres).
// Certified bracket for comparison: max(uA, uB) <= u <= min(uA + uB, 1).
// ---------------------------------------------------------------------------
float ptUnionCoupled( float uA, float uB, float uBatA, float uAatB )
{
    float den = 1.0 - uBatA * uAatB;
    float al = (1.0 - uBatA) / den, be = (1.0 - uAatB) / den;
    return clamp( al * uA + be * uB, max( uA, uB ), 1.0 );
}

#endif // SDF_GEMS_HAS_PT_HARMONIC
