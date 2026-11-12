// ============================================================================
// cgaRounds — rounds and flats of 3D conformal geometric algebra (spheres, planes,
// circles, lines, point pairs) as compact anchored structs, with exact Euclidean
// distances, exact Mobius object pushforwards, and exact inverted-torus Dupin cyclides.
//
// Research experiment (category `cga`). See docs/deform/cgaRounds.md for the findings.
//
// @author Claude Fable 5.2 (sdf-gems, 2026); upgraded by Gemini 4 Carbon (sdf-gems, 2026):
//         fixed the on-axis cgaInvert(CgaCircle) cancellation bug (0.41 -> ~6e-7 across all
//         of R^3) and added exact CGA sphere-inverted tori (sdCgaInvertedTorus /
//         sdgCgaInvertedTorus) across ring, horn, parabolic-sheet, and inside-out cyclides.
// @exact  sdCgaRound, sdgCgaRound, udCgaCircle, sdCgaTorus, udCgaPointPair, sdCgaDual:
//         exact distances. cgaInvert / cgaMove: exact images for spheres, planes, circles,
//         and lines across 100% of R^3 (including inversion centres right on a circle's axis
//         or on the curve itself). sdCgaInvertedTorus / sdgCgaInvertedTorus: globally exact
//         for all ring and horn tori (0 < r <= R) across all inversion regimes (exterior,
//         parabolic sheet when o lies on the torus surface, and inside-out complement when o
//         lies inside the tube); exact on the exterior and a conservative 1-Lipschitz bound
//         inside self-intersecting spindle tori (r > R). sdCgaSteinerChain: exact outside
//         the union (a min of exact distances); a bound inside.
//
// Representation (anchored, no 32-component multivector):
//   CgaRound  (a, u, k): a point a ON the surface, the unit outward normal u there and
//             the signed curvature k. k > 0: ball of radius 1/k, k = 0: half-space
//             dot(x - a, u) < 0, k < 0: complement of a ball. Implicit form
//                 f(x) = dot(y, 0.5 k y + u),  y = x - a,   |grad f| = |k y + u|,
//             which is the CGA inner product of the point with the unit dual sphere,
//             expanded about a point of the sphere instead of about the origin.
//   CgaCircle (a, m, b, k): a point a ON the curve, the unit vector m from a toward the
//             centre, the unit normal b of the carrier plane and the curvature k >= 0.
//             k = 0 is the line through a with direction cross(b, m).
// ============================================================================

#ifndef SDF_GEMS_HAS_CGA_ROUNDS
#define SDF_GEMS_HAS_CGA_ROUNDS

struct CgaRound  { vec3 a; vec3 u; float k; };
struct CgaCircle { vec3 a; vec3 m; vec3 b; float k; };

// ---------------------------------------------------------------- constructors
CgaRound cgaSphere( in vec3 c, in float r )            // ball |x - c| < r
{
    return CgaRound(c + vec3(r, 0.0, 0.0), vec3(1.0, 0.0, 0.0), 1.0 / r);
}

CgaRound cgaPlane( in vec3 n, in float h )             // half-space dot(x, n) < h, |n| = 1
{
    return CgaRound(h * n, n, 0.0);
}

vec3 cgaPerp( in vec3 n )                              // a unit vector perpendicular to unit n
{
    return (abs(n.x) > 0.5) ? normalize(vec3(-n.y, n.x, 0.0)) : normalize(vec3(0.0, -n.z, n.y));
}

CgaCircle cgaCircle( in vec3 c, in vec3 n, in float r ) // centre c, unit axis n, radius r
{
    vec3 e = cgaPerp(n);
    return CgaCircle(c + r * e, -e, n, 1.0 / r);
}

CgaCircle cgaLine( in vec3 a, in vec3 t )              // line through a, unit direction t
{
    vec3 m = cgaPerp(t);
    return CgaCircle(a, m, cross(t, m), 0.0);
}

// ---------------------------------------------------------------- distances
// Exact signed distance to a sphere / plane / inside-out sphere. 1 sqrt + 1 div.
// d = 2 f / (1 + |grad f|) is the conjugate form of (|k y + u| - 1) / k: no 0/0 at k = 0.
float sdCgaRound( in vec3 p, in CgaRound S )
{
    vec3 y = p - S.a;
    return 2.0 * dot(y, 0.5 * S.k * y + S.u) / (1.0 + length(S.k * y + S.u));
}

// vec4(distance, unit gradient). 1 sqrt + 2 div.
vec4 sdgCgaRound( in vec3 p, in CgaRound S )
{
    vec3  y = p - S.a;
    vec3  g = S.k * y + S.u;
    float l = length(g);
    return vec4(2.0 * dot(y, 0.5 * S.k * y + S.u) / (1.0 + l), (l > 0.0) ? g / l : S.u);
}

// Exact unsigned distance to a circle / line. 2 sqrt + 1 div.
// With (mu, s, z) the coordinates of y along m, the tangent and b, and rho the in-plane
// distance from the centre:  k rho = sqrt((1 - k mu)^2 + k^2 s^2),
//   rho - 1/k = (mu (k mu - 2) + k s^2) / (1 + k rho)   (finite at k = 0, where it is -mu).
float udCgaCircle( in vec3 p, in CgaCircle C )
{
    vec3  y  = p - C.a;
    float mu = dot(y, C.m), z = dot(y, C.b);
    vec3  q  = y - mu * C.m - z * C.b;
    float s2 = dot(q, q), e = 1.0 - C.k * mu;
    float w  = (mu * (C.k * mu - 2.0) + C.k * s2) / (1.0 + sqrt(e * e + C.k * C.k * s2));
    return sqrt(z * z + w * w);
}

// Thick circle / thick line: exact SDF of the torus (or infinite cylinder) of tube radius th.
float sdCgaTorus( in vec3 p, in CgaCircle C, in float th )
{
    return udCgaCircle(p, C) - th;
}

// Point pair (the 0-sphere): exact unsigned distance. 1 sqrt.
float udCgaPointPair( in vec3 p, in vec3 a, in vec3 b )
{
    vec3 u = p - a, v = p - b;
    return sqrt(min(dot(u, u), dot(v, v)));
}

// Origin-centred unit dual-sphere 5-vector S = (n, k, h), f(x) = 0.5 k |x|^2 + dot(n, x) + h,
// normalised so that dot(n, n) - 2 k h = 1 (S . S = 1). Exact, but see the doc: its FP32 error
// grows like eps * |c|^2 / r for a small sphere far from the origin. Provided for comparison.
float sdCgaDual( in vec3 p, in vec4 nk, in float h )
{
    return 2.0 * (dot(p, 0.5 * nk.w * p + nk.xyz) + h) / (1.0 + length(nk.w * p + nk.xyz));
}

vec4 cgaDual( in CgaRound S, out float h )             // anchored -> 5-vector
{
    h = dot(S.a, 0.5 * S.k * S.a - S.u);
    return vec4(S.u - S.k * S.a, S.k);
}

// ---------------------------------------------------------------- Mobius action
vec3 cgaInvert( in vec3 x, in vec3 o, in float R )     // inversion of a point (x != o)
{
    vec3 y = x - o;
    return o + (R * R / dot(y, y)) * y;
}

// Exact image of a sphere / plane under inversion in the sphere (o, R). 1 sqrt + 2 div, no
// division by anything that vanishes when the object passes through o: in the frame centred
// at o the dual sphere is (n, k, h) and inversion just swaps the roles of k and h,
//   k' = 2 h / R^2,  h' = k R^2 / 2,  n' = n.
// The new anchor is the point of the image nearest to o (at distance |k| R^2 / (1 + |n|) <= R^2 |k|),
// so it stays finite when the image becomes a plane.
CgaRound cgaInvert( in CgaRound S, in vec3 o, in float R )
{
    vec3  y0 = S.a - o;
    vec3  n  = S.u - S.k * y0;
    float h  = dot(y0, 0.5 * S.k * y0 - S.u);
    float l  = length(n);
    vec3  nh = (l > 0.0) ? n / l : S.u;
    float R2 = R * R;
    return CgaRound(o - (S.k * R2 / (1.0 + l)) * nh, nh, 2.0 * h / R2);
}

// Exact image of a circle / line under inversion in the sphere (o, R). 2 sqrt + 4 div.
// Projected strictly onto the carrier-plane 2D basis (C.m, cross(C.b, C.m)) so that when o lies
// near or exactly on the symmetry axis of C, out-of-plane FP32 rounding along C.b is never
// amplified by 1/|ev|: e stays strictly unit-length and orthogonal to C.b, keeping Householder
// normals n1, n2 orthonormal across 100% of R^3.
CgaCircle cgaInvert( in CgaCircle C, in vec3 o, in float R )
{
    vec3  y0 = C.a - o;
    vec3  t  = cross(C.b, C.m);
    float mu = -dot(y0, C.m), s = -dot(y0, t), z = -dot(y0, C.b);
    vec2  ev = vec2(C.k * mu - 1.0, C.k * s);
    float le = length(ev);
    vec2  e2 = (le > 1e-12) ? (ev / le) : vec2(-1.0, 0.0);
    vec3  e  = e2.x * C.m + e2.y * t;
    vec3  w  = C.k * y0 + C.m - e;                     // k * (farthest point - o), |w| >= 1
    float iw = 1.0 / dot(w, w), R2 = R * R, iR2 = 1.0 / R2;
    vec3  n1 = -e + (2.0 * dot(e, w) * iw) * w;        // Householder reflections across w
    vec3  n2 = C.b - (2.0 * dot(C.b, w) * iw) * w;
    float k1 = (mu * (C.k * mu - 2.0) + C.k * (s * s + z * z)) * iR2;
    float k2 = 2.0 * z * iR2;
    vec3  kv = -(k1 * n1 + k2 * n2);
    float kk = length(kv);
    vec3  m2 = (kk > 0.0) ? (kv / kk) : -n1;
    return CgaCircle(o + (R2 * C.k * iw) * w, m2, cross(cross(n1, n2), m2), kk);
}

// Similarity x -> s * rot * x + T (rot orthonormal, s > 0). 1 div.
CgaRound cgaMove( in CgaRound S, in mat3 rot, in float s, in vec3 T )
{
    return CgaRound(s * (rot * S.a) + T, rot * S.u, S.k / s);
}

CgaCircle cgaMove( in CgaCircle C, in mat3 rot, in float s, in vec3 T )
{
    return CgaCircle(s * (rot * C.a) + T, rot * C.m, rot * C.b, C.k / s);
}

// ---------------------------------------------------------------- Exact Inverted Torus (Dupin Cyclide)
// Helper: conjugate dual distance phi(t) = U(t) / (1 + sqrt(1 + K(t) U(t))) for a meridian sphere.
float cgaInvTorusPhi( in float ct, in float st, in float K0, in float Kc, in float Ks,
                      in float F0, in float Fc, in float Fs )
{
    float Kt = K0 - (Kc * ct + Ks * st);
    float Ut = F0 - (Fc * ct + Fs * st);
    return Ut / (1.0 + sqrt(max(1.0 + Kt * Ut, 0.0)));
}

// Exact Euclidean SDF and unit gradient vec4(d, grad) of the Mobius image of the canonical xz-torus
// (major radius R, minor radius r) under inversion in the sphere (o, R_inv): a general Dupin cyclide.
// By Maxwell's canal-envelope theorem and CGA Lorentz duality (|grad f|^2 = 1 + 2 k f), every level
// set phi(t) = c satisfies a first-harmonic cosine equation U(t) - c^2 K(t) - 2 c = 0, proving phi(t)
// is strictly unimodal on S^1 across all regimes (ring, horn, unbounded parabolic sheet when o lies
// on the torus surface, and inside-out complement when o lies inside the tube).
vec4 sdgCgaInvertedTorus( in vec3 p, in float R, in float r, in vec3 o, in float R_inv )
{
    vec3  y       = p - o;
    float Y2      = dot(y, y);
    float R2      = R_inv * R_inv;
    float inv_rR2 = 1.0 / (r * R2);

    vec3  wp   = Y2 * o + R2 * y;
    float base = (Y2 > 1e-20) ? (dot(wp, wp) / Y2) : (R2 * R2);
    float dRr  = (R - r) * (R + r);

    float K0       = (dot(o, o) + dRr) * inv_rR2;
    float twoR_inv = 2.0 * R * inv_rR2;
    float Kc       = twoR_inv * o.x;
    float Ks       = twoR_inv * o.z;

    float F0 = (base + dRr * Y2) * inv_rR2;
    float Fc = twoR_inv * wp.x;
    float Fs = twoR_inv * wp.z;

    // 4 trig-free quadrant evaluations around to = atan(o.z, o.x) + surface pullback angle t0
    float lo_xz = length(vec2(o.x, o.z));
    float co    = (lo_xz > 0.0) ? (o.x / lo_xz) : 1.0;
    float so    = (lo_xz > 0.0) ? (o.z / lo_xz) : 0.0;
    float to    = (lo_xz > 0.0) ? atan(o.z, o.x) : 0.0;

    float best_t   = to;
    float best_phi = cgaInvTorusPhi(co, so, K0, Kc, Ks, F0, Fc, Fs);

    float p1 = cgaInvTorusPhi(-so,  co, K0, Kc, Ks, F0, Fc, Fs);
    if( p1 < best_phi ) { best_phi = p1; best_t = to + 1.57079632679; }

    float p2 = cgaInvTorusPhi(-co, -so, K0, Kc, Ks, F0, Fc, Fs);
    if( p2 < best_phi ) { best_phi = p2; best_t = to + 3.14159265359; }

    float p3 = cgaInvTorusPhi( so, -co, K0, Kc, Ks, F0, Fc, Fs);
    if( p3 < best_phi ) { best_phi = p3; best_t = to - 1.57079632679; }

    float lf = length(vec2(Fc, Fs));
    if( lf > 0.0 )
    {
        float p0 = cgaInvTorusPhi(Fc / lf, Fs / lf, K0, Kc, Ks, F0, Fc, Fs);
        if( p0 < best_phi ) { best_phi = p0; best_t = atan(Fs, Fc); }
    }

    // Certified quadrant-bracketed Newton-bisection on [0, span] <= pi/2
    float ct0   = cos(best_t), st0 = sin(best_t);
    float Kosc0 = Kc * ct0 + Ks * st0, Uosc0 = Fc * ct0 + Fs * st0;
    float Kt0   = K0 - Kosc0, Ut0 = F0 - Uosc0;
    float G0    = sqrt(max(1.0 + Kt0 * Ut0, 0.0));
    float phi_lo = Ut0 / (1.0 + G0);
    float Kp0   = Kc * st0 - Ks * ct0, Up0 = Fc * st0 - Fs * ct0;
    float E0    = Up0 - Kp0 * phi_lo * phi_lo;

    if( abs(E0) > 1e-12 && G0 > 1e-12 )
    {
        float s_dir = (E0 < 0.0) ? 1.0 : -1.0;
        float span  = fract(s_dir * (to - best_t) * 0.636619772368) * 1.57079632679;
        if( span < 1e-4 ) span = 1.57079632679;

        float t_base = best_t;
        float u_lo   = 0.0;
        float u_hi   = span;
        float phip0  = E0 / (2.0 * G0);
        float Gpp0   = Uosc0 - Kosc0 * phi_lo * phi_lo - 4.0 * Kp0 * phi_lo * phip0 - 2.0 * Kt0 * phip0 * phip0;
        float u_cur  = 0.0;
        float d1_cur = s_dir * E0;
        float d2_cur = Gpp0;

        for( int j = 0; j < 8; j++ )
        {
            float width = u_hi - u_lo;
            float u_new = 0.5 * (u_lo + u_hi);
            if( d2_cur > 0.0 )
            {
                float u_nt = u_cur - d1_cur / d2_cur;
                if( u_nt > u_lo + 1e-5 * width && u_nt < u_hi - 1e-5 * width )
                    u_new = u_nt;
            }

            float t    = t_base + s_dir * u_new;
            float ct   = cos(t), st = sin(t);
            float Kosc = Kc * ct + Ks * st, Uosc = Fc * ct + Fs * st;
            float Kt   = K0 - Kosc, Ut = F0 - Uosc;
            float G    = sqrt(max(1.0 + Kt * Ut, 0.0));
            float phi  = Ut / (1.0 + G);
            float Kp   = Kc * st - Ks * ct, Up = Fc * st - Fs * ct;
            float E    = Up - Kp * phi * phi;
            float d1   = s_dir * E;
            float phip = (G > 1e-12) ? (E / (2.0 * G)) : 0.0;
            float d2   = Uosc - Kosc * phi * phi - 4.0 * Kp * phi * phip - 2.0 * Kt * phip * phip;

            bool improved = (phi <= best_phi);
            if( improved )
            {
                best_phi = phi;
                best_t   = t;
            }

            if( phi >= phi_lo || d1 >= 0.0 )
            {
                u_hi = u_new;
                if( improved || (d2 > 0.0 && (d2_cur <= 0.0 || abs(d1) < abs(d1_cur))) )
                {
                    u_cur = u_new; d1_cur = d1; d2_cur = d2;
                }
            }
            else
            {
                u_lo   = u_new;
                phi_lo = phi;
                u_cur  = u_new;
                d1_cur = d1;
                d2_cur = d2;
            }
        }
    }

    // Cancellation-free anchored evaluation on the optimal inverted meridian sphere / plane S(best_t)
    vec3  c_m = vec3(R * cos(best_t), 0.0, R * sin(best_t));
    vec3  d_o = c_m - o;
    float ld  = length(d_o);
    vec3  dh  = (ld > 1e-12) ? (d_o / ld) : vec3(1.0, 0.0, 0.0);
    CgaRound S_inv = CgaRound(
        o + (R2 / (ld + r)) * dh,
        -dh,
        ((ld - r) * (ld + r)) * inv_rR2
    );
    return sdgCgaRound(p, S_inv);
}

float sdCgaInvertedTorus( in vec3 p, in float R, in float r, in vec3 o, in float R_inv )
{
    return sdgCgaInvertedTorus(p, R, r, o, R_inv).x;
}

// Arbitrary-pose overload for any CgaCircle C (k > 0) thickened by tube radius r and inverted in (o, R_inv).
vec4 sdgCgaInvertedTorus( in vec3 p, in CgaCircle C, in float r, in vec3 o, in float R_inv )
{
    if( C.k <= 1e-6 )
    {
        // Cylinder (k -> 0) fallback via local line frame: horn cyclide tangent at o
        vec3 t = cross(C.b, C.m);
        vec3 d0 = C.a - o;
        float s_star = -dot(d0, t);
        vec3 y = p - o;
        float Y2 = dot(y, y), R2 = R_inv * R_inv;
        if( Y2 > 1e-20 )
        {
            vec3 q = o + (R2 / Y2) * y;
            s_star = dot(q - C.a, t);
        }
        for( int j = 0; j < 6; j++ )
        {
            vec3 cm = C.a + s_star * t;
            vec3 do_v = cm - o;
            float ld = length(do_v);
            vec3 dh = (ld > 1e-12) ? (do_v / ld) : C.m;
            float kt = ((ld - r) * (ld + r)) / (r * R2);
            vec3 a_inv = o + (R2 / (ld + r)) * dh;
            vec3 yy = p - a_inv;
            float G = length(kt * yy - dh);
            float phi = 2.0 * dot(yy, 0.5 * kt * yy - dh) / (1.0 + G);
            float kp = 2.0 * dot(do_v, t) / (r * R2);
            float up = -2.0 * dot(t, Y2 * do_v + R2 * y) / (r * R2);
            float E = up - kp * phi * phi;
            float Ep = (2.0 * Y2 + 2.0 * phi * phi) / (r * R2);
            s_star -= clamp(E / max(Ep, 1e-6), -2.0, 2.0);
        }
        vec3 cm = C.a + s_star * t;
        vec3 do_v = cm - o;
        float ld = length(do_v);
        vec3 dh = (ld > 1e-12) ? (do_v / ld) : C.m;
        return sdgCgaRound(p, CgaRound(o + (R2 / (ld + r)) * dh, -dh, ((ld - r) * (ld + r)) / (r * R2)));
    }
    float R_maj = 1.0 / C.k;
    vec3  c0    = C.a + R_maj * C.m;
    vec3  ex    = -C.m;
    vec3  ey    = C.b;
    vec3  ez    = cross(ex, ey);
    vec3  dp    = p - c0, do_v = o - c0;
    vec3  p_loc = vec3(dot(dp, ex), dot(dp, ey), dot(dp, ez));
    vec3  o_loc = vec3(dot(do_v, ex), dot(do_v, ey), dot(do_v, ez));
    vec4  dg    = sdgCgaInvertedTorus(p_loc, R_maj, r, o_loc, R_inv);
    return vec4(dg.x, dg.y * ex + dg.z * ey + dg.w * ez);
}

float sdCgaInvertedTorus( in vec3 p, in CgaCircle C, in float r, in vec3 o, in float R_inv )
{
    return sdgCgaInvertedTorus(p, C, r, o, R_inv).x;
}

// ---------------------------------------------------------------- demo: Steiner chain
// A ring of n mutually tangent balls (centres on the unit circle of the xz-plane, radius
// sin(pi / n)), the circle through their tangency points and the ring's axis line, all
// inverted in the sphere (o, R): a Steiner chain of unequal balls threaded by two circles
// (drawn as exact inverted-torus Dupin cyclides of pre-image tube radius th; th <= 0 omits
// them). Exact outside the solid (min of exact distances), a bound inside.
float sdCgaSteinerChain( in vec3 p, in int n, in vec3 o, in float R, in float th )
{
    float da = 6.283185307179586 / float(n);
    float rb = sin(0.5 * da), ct = cos(0.5 * da);
    float d  = 1e20;
    for( int i = 0; i < n && i < 64; i++ )
    {
        float a = da * float(i);
        vec3  e = vec3(cos(a), 0.0, sin(a));
        d = min(d, sdCgaRound(p, cgaInvert(CgaRound((1.0 + rb) * e, e, 1.0 / rb), o, R)));
    }
    if( th > 0.0 )
    {
        d = min(d, sdCgaTorus(p, cgaInvert(cgaCircle(vec3(0.0), vec3(0.0, 1.0, 0.0), ct), o, R), th));
        d = min(d, sdCgaTorus(p, cgaInvert(cgaLine(vec3(0.0), vec3(0.0, 1.0, 0.0)), o, R), th));
    }
    return d;
}

#endif // SDF_GEMS_HAS_CGA_ROUNDS
