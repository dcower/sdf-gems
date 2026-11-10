// ============================================================================
// cgaRounds — rounds and flats of 3D conformal geometric algebra (spheres, planes,
// circles, lines, point pairs) as compact anchored structs, with exact Euclidean
// distances and the exact Mobius action on the OBJECT (inversion + similarities).
//
// Research experiment (category `cga`). See docs/cga/cgaRounds.md for the findings.
//
// @author Claude Fable 5.2 (sdf-gems, 2026). The geometry is classical (inversive /
//         conformal geometric algebra); see the doc for prior art.
// @exact  sdCgaRound, sdgCgaRound, udCgaCircle, sdCgaTorus, udCgaPointPair, sdCgaDual:
//         exact distances. cgaInvert / cgaMove: exact images (zero-thickness objects,
//         balls and half-spaces only). sdCgaSteinerChain: exact outside the union
//         (a min of exact distances); a bound inside.
//         MEASURED FAILURE: cgaInvert(CgaCircle) is ill-conditioned when o is near the
//         circle's axis (error ~ 1e-7 / distance of o from the axis; wrong ON the axis).
//         NOT covered: the Mobius image of a thick circle (torus) is a Dupin cyclide,
//         not sdCgaTorus of the image circle.
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

// Exact image of a circle / line under inversion in the sphere (o, R). 3 sqrt + 4 div.
// CAVEAT (measured): NOT reliable when o is on or near the circle's axis, see the header.
// The image is anchored at the image of the point of C FARTHEST from o (the origin o itself
// for a line), so the anchor stays within R^2 k of o. The image is the intersection of the
// images of the carrier plane and of the sphere through C centred in that plane; both pass
// through the new anchor with orthogonal normals n1, n2 and curvatures k1, k2, and the
// curvature vector of their intersection curve is -(k1 n1 + k2 n2).
CgaCircle cgaInvert( in CgaCircle C, in vec3 o, in float R )
{
    vec3  y0 = C.a - o;
    float mu = -dot(y0, C.m), z = -dot(y0, C.b);
    vec3  q  = -y0 - mu * C.m - z * C.b;
    vec3  ev = (C.k * mu - 1.0) * C.m + C.k * q;
    float le = length(ev);
    vec3  e  = (le > 0.0) ? ev / le : -C.m;
    vec3  w  = C.k * y0 + C.m - e;                     // k * (farthest point - o), |w| >= 1
    float iw = 1.0 / dot(w, w), R2 = R * R, iR2 = 1.0 / R2;
    vec3  n1 = -e + (2.0 * dot(e, w) * iw) * w;        // Householder reflections across w
    vec3  n2 = C.b - (2.0 * dot(C.b, w) * iw) * w;
    float k1 = (C.k * dot(y0, y0) + 2.0 * dot(C.m, y0)) * iR2;
    float k2 = -2.0 * dot(C.b, y0) * iR2;
    vec3  kv = -(k1 * n1 + k2 * n2);
    float kk = length(kv);
    vec3  m2 = (kk > 0.0) ? kv / kk : -n1;
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

// ---------------------------------------------------------------- demo: Steiner chain
// A ring of n mutually tangent balls (centres on the unit circle of the xz-plane, radius
// sin(pi / n)), the circle through their tangency points and the ring's axis line, all
// inverted in the sphere (o, R): a Steiner chain of unequal balls threaded by two circles
// (drawn as tori of tube radius th; th <= 0 omits them). Exact outside the solid (min of
// exact distances), a bound inside. In production hoist the n object transforms out of the
// per-pixel code (uniforms / vertex stage); they do not depend on p.
float sdCgaSteinerChain( in vec3 p, in int n, in vec3 o, in float R, in float th )
{
    float da = 6.283185307179586 / float(n);
    float rb = sin(0.5 * da), ct = cos(0.5 * da);
    float d  = 1e20;
    for( int i = 0; i < n && i < 64; i++ )
    {
        float a = da * float(i);
        vec3  e = vec3(cos(a), 0.0, sin(a));
        // anchor each ball at its outermost point (any point of the sphere would do)
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
