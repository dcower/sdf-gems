// ============================================================================
// msSteinerTube — Closed-form intrinsic volumes (V, A, M, chi) of SDF primitives,
//                 Steiner / Weyl tube polynomials, offset bodies, shells and
//                 closed-form "glow" integrals of radial kernels g(d(x)).
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @source Steiner (1840), Weyl (1939), Federer (1959): textbook integral geometry;
//         packaged here for an SDF library and verified against GPU Monte Carlo
//         histograms of the literal repo SDFs (see docs/duality/msSteinerTube.md).
// @exact  The measures of each primitive are closed-form (verified to Monte Carlo
//         accuracy, ~1e-4 relative). The tube polynomial is an identity ONLY for
//         -reach_in < r < reach_out of the body (all r >= 0 for convex bodies;
//         r <= 0 fails immediately for bodies with sharp convex edges).
//         Outside that range it is NOT a bound in either direction (measured: it
//         under-estimates by 1.7% for a torus whose hole has closed and over-estimates
//         by 27% for two linked tori whose tubes overlap).
// ============================================================================
//
// Convention: every ms* constructor returns  m = vec4(V, A, M, chi)  with
//   V   = volume,
//   A   = surface area,
//   M   = integral of the mean curvature H = (k1 + k2) / 2 over the surface
//         (a sharp convex edge of length L and exterior dihedral angle t adds t * L / 2;
//          for convex bodies M = 2 pi * mean width),
//   chi = Euler characteristic of the solid (ball 1, solid torus 0, spherical shell 2)
//         = (1 / 4 pi) * integral of the Gauss curvature.
// Tube (Steiner / Weyl) formula for the exact SDF d of the body:
//   vol { d < r } = V + A r + M r^2 + (4 pi / 3) chi r^3
//   area{ d = r } = A + 2 M r + 4 pi chi r^2

#ifndef SDF_GEMS_MS_STEINER_TUBE_GLSL
#define SDF_GEMS_MS_STEINER_TUBE_GLSL

// ---------------------------------------------------------------------------
// 1. Tube polynomial, offsets, shells, unions
// ---------------------------------------------------------------------------

// Measures of the offset body { d < r } (Taylor shift of the Steiner polynomial).
// Valid for -reach_in < r < reach_out.
vec4 msOffset( vec4 m, float r )
{
    const float PI4 = 12.566370614359172;
    float k = PI4 * m.w;
    return vec4(m.x + r * (m.y + r * (m.z + r * k / 3.0)),
                m.y + r * (2.0 * m.z + r * k),
                m.z + r * k,
                m.w);
}

// vol{ d < r }
float msTubeVolume( vec4 m, float r )
{
    return m.x + r * (m.y + r * (m.z + r * 4.1887902047863905 * m.w));
}

// area{ d = r }  ( = d/dr msTubeVolume )
float msTubeArea( vec4 m, float r )
{
    return m.y + r * (2.0 * m.z + r * 12.566370614359172 * m.w);
}

// Shell { |d| < t } (the "onion" operator abs(d) - t). Needs t < min(reach_in, reach_out).
// The inner wall is counted with reversed orientation: A adds, M subtracts, chi doubles.
vec4 msOnion( vec4 m, float t )
{
    vec4 o = msOffset(m, t), i = msOffset(m, -t);
    return vec4(o.x - i.x, o.y + i.y, o.z - i.z, o.w + i.w);
}

// Union of two bodies whose closures are disjoint. The tube formula of the union holds
// for r < min(reach_a, reach_b, gap / 2).
vec4 msUnionDisjoint( vec4 a, vec4 b )
{
    return a + b;
}

// ---------------------------------------------------------------------------
// 2. Primitives (same parameters as the matching sd* functions in primitives3d/)
// ---------------------------------------------------------------------------

// sdSphere(p, s)
vec4 msSphere( float s )
{
    const float PI4 = 12.566370614359172;
    return vec4(PI4 * s * s * s / 3.0, PI4 * s * s, PI4 * s, 1.0);
}

// sdBox(p, b), b = half extents.  M = pi * (sum of the three full edge lengths).
vec4 msBox( vec3 b )
{
    return vec4(8.0 * b.x * b.y * b.z,
                8.0 * (b.x * b.y + b.y * b.z + b.z * b.x),
                6.283185307179586 * (b.x + b.y + b.z),
                1.0);
}

// sdRoundBox(p, b, r): box of half extents b - r, offset by r.
vec4 msRoundBox( vec3 b, float r )
{
    return msOffset(msBox(b - r), r);
}

// Capsule: segment of length len swept by a ball of radius r.
vec4 msCapsule( float len, float r )
{
    const float PI = 3.141592653589793;
    return vec4(PI * r * r * (len + r * 4.0 / 3.0),
                2.0 * PI * r * (len + 2.0 * r),
                PI * (len + 4.0 * r),
                1.0);
}

// sdCapsule(p, a, b, r)
vec4 msCapsule( vec3 a, vec3 b, float r )
{
    return msCapsule(length(b - a), r);
}

// sdCappedCylinder(p, h, r): half height h, radius r.  M = pi * (2h) + pi^2 r.
vec4 msCylinder( float h, float r )
{
    const float PI = 3.141592653589793;
    return vec4(2.0 * PI * r * r * h,
                2.0 * PI * r * (2.0 * h + r),
                PI * (2.0 * h + PI * r),
                1.0);
}

// sdRoundedCylinder(p, ra, rb, h): cylinder of radius 2 ra - rb, half height h, offset by rb.
vec4 msRoundedCylinder( float ra, float rb, float h )
{
    return msOffset(msCylinder(h, 2.0 * ra - rb), rb);
}

// sdCappedCone(p, h, r1, r2): half height h, radius r1 at y = -h and r2 at y = +h.
// Lateral surface: integral of H = pi * height. Rims: exterior angles pi/2 +- alpha.
vec4 msCappedCone( float h, float r1, float r2 )
{
    const float PI = 3.141592653589793;
    float hh = 2.0 * h, dr = r1 - r2, sl = sqrt(hh * hh + dr * dr), al = atan(dr, hh);
    return vec4(PI * hh * (r1 * r1 + r1 * r2 + r2 * r2) / 3.0,
                PI * (r1 * r1 + r2 * r2 + (r1 + r2) * sl),
                PI * (hh + 0.5 * PI * (r1 + r2) + al * dr),
                1.0);
}

// sdCone(p, c, h): c = (sin, cos) of the half apex angle, height h.
vec4 msCone( vec2 c, float h )
{
    return msCappedCone(0.5 * h, h * c.x / c.y, 0.0);
}

// sdTorus(p, t): t = (major radius, minor radius). Tube formula valid for -t.y < r < t.x - t.y.
vec4 msTorus( vec2 t )
{
    const float PI2SQ = 19.739208802178716; // 2 pi^2
    return vec4(PI2SQ * t.x * t.y * t.y, 2.0 * PI2SQ * t.x * t.y, PI2SQ * t.x, 0.0);
}

// Solid tube of radius a around ANY closed space curve of length len (Weyl: independent of the
// shape of the curve, knotted or not), valid while a + r is below the reach of the curve.
vec4 msTubeCurve( float len, float a )
{
    const float PI = 3.141592653589793;
    return vec4(PI * a * a * len, 2.0 * PI * a * len, PI * len, 0.0);
}

// sdLink(p, le, r1, r2): stadium curve of length 2 pi r1 + 4 le, tube radius r2.
vec4 msLink( float le, float r1, float r2 )
{
    return msTubeCurve(6.283185307179586 * r1 + 4.0 * le, r2);
}

// Right prism of height hgt over ANY convex polygon / convex region of area f and perimeter per.
// (vertical edges: total exterior angle 2 pi -> pi * hgt;  rims: (pi / 4) * per each.)
vec4 msPrism( float f, float per, float hgt )
{
    const float PI = 3.141592653589793;
    return vec4(f * hgt, 2.0 * f + per * hgt, PI * (hgt + 0.5 * per), 1.0);
}

// sdHexPrism(p, h): apothem h.x, half height h.y.
vec4 msHexPrism( vec2 h )
{
    const float S3 = 1.7320508075688772;
    return msPrism(2.0 * S3 * h.x * h.x, 4.0 * S3 * h.x, 2.0 * h.y);
}

// sdSpheroid(p, ab): equatorial radius ab.x (xz plane), polar semi-axis ab.y.
vec4 msSpheroid( vec2 ab )
{
    const float PI = 3.141592653589793;
    float a = ab.x, c = ab.y, a2 = a * a, c2 = c * c;
    float s = (c2 - a2) / a2;                 // > 0 prolate, < 0 oblate
    float fa, fm;                             // A = 2 pi a^2 (1 + fa),  M = 2 pi (c + a fm)
    if( abs(s) < 1e-3 )
    {
        fa = 1.0 + s * (2.0 / 3.0 - s * (2.0 / 15.0));           // series about the sphere
        fm = 1.0 - s * (1.0 / 6.0 - s * 0.075);
    }
    else if( s > 0.0 )
    {
        float e = sqrt(1.0 - a2 / c2), x = sqrt(s);
        fa = c * asin(e) / (a * e);
        fm = asinh(x) / x;
    }
    else
    {
        float e = sqrt(1.0 - c2 / a2);
        fa = (c2 / a2) * atanh(e) / e;
        fm = asin(e) / e;
    }
    return vec4(4.0 * PI * a2 * c / 3.0, 2.0 * PI * a2 * (1.0 + fa), 2.0 * PI * (c + a * fm), 1.0);
}

// ---------------------------------------------------------------------------
// 3. Closed-form integrals of a kernel g(d(x)) over the EXTERIOR { d > 0 }:
//      I = int_0^inf g(r) (A + 2 M r + 4 pi chi r^2) dr.
//    Identities for convex bodies (chi = 1). For non-convex bodies they hold only when the
//    kernel is supported inside the reach (msGlowTent with R < reach); otherwise they are
//    approximations of either sign, with an error that vanishes as sigma / reach -> 0.
// ---------------------------------------------------------------------------

// g = exp(-d / sigma)
float msGlowExp( vec4 m, float sigma )
{
    return sigma * (m.y + sigma * (2.0 * m.z + sigma * 25.132741228718345 * m.w));
}

// g = exp(-d^2 / (2 sigma^2))
float msGlowGauss( vec4 m, float sigma )
{
    const float K = 1.2533141373155003; // sqrt(pi / 2)
    return sigma * (K * m.y + sigma * (2.0 * m.z + sigma * K * 12.566370614359172 * m.w));
}

// g = 1 / (1 + (d / sigma)^2)^2   (squared Lorentzian; a plain 1 / (1 + d^2) diverges in 3D)
float msGlowInvSq( vec4 m, float sigma )
{
    return sigma * (0.7853981633974483 * m.y + sigma * (m.z + sigma * 9.869604401089358 * m.w));
}

// g = (1 + d / sigma)^(-n),  n > 3
float msGlowPow( vec4 m, float sigma, float n )
{
    float a = 1.0 / (n - 1.0), b = a / (n - 2.0), c = 2.0 * b / (n - 3.0);
    return sigma * (a * m.y + sigma * (2.0 * b * m.z + sigma * c * 12.566370614359172 * m.w));
}

// g = max(1 - d / R, 0)   (compact support: exact for any body with reach_out >= R)
float msGlowTent( vec4 m, float R )
{
    return R * (0.5 * m.y + R * (m.z / 3.0 + R * 1.0471975511965976 * m.w));
}

// ---------------------------------------------------------------------------
// 4. Convex-body corollaries
// ---------------------------------------------------------------------------

// Mean width (average over all directions of the extent of the body), convex bodies: M / (2 pi).
// Additive under Minkowski sums.
float msMeanWidth( vec4 m )
{
    return m.z * 0.15915494309189535;
}

// Cauchy: the mean orthographic shadow area of a convex body over all directions is A / 4.
float msMeanProjectedArea( vec4 m )
{
    return 0.25 * m.y;
}

// Mean chord length of a convex body under uniform random lines: 4 V / A.
float msMeanChord( vec4 m )
{
    return 4.0 * m.x / m.y;
}

#endif // SDF_GEMS_MS_STEINER_TUBE_GLSL
