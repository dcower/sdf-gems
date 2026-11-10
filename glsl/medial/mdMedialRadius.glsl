// ============================================================================
// mdMedialRadius — Closed-Form Medial Radius Fields paired with exact SDFs
//   mdSphere, mdBox, mdRoundBox, mdCapsule, mdCylinder, mdTorus, mdTorusExterior,
//   mdSpheroid, mdCone, mdEllipsoidFoot, mdEllipsoidFromSdg
//   + operators mdOffset, mdOffsetIsExact, mdThickness, mdSkeleton
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  d: Yes (same exact SDFs as the sd* primitives, re-derived here so the file is self-contained).
//         rho: closed form, verified against two independent ground truths (see docs/medial/mdMedialRadius.md).
//         mdSkeleton: NOT a distance bound from below (it is an upper bound on the distance to the medial axis).
//
// Every md* function returns vec2(d, rho):
//   d   = exact signed distance to the surface (negative inside),
//   rho = medial radius of the foot point x of p: the radius of the largest ball inside the body that
//         touches the surface at x (the distance from x, along the inward normal, to the medial axis).
//         rho is constant along each normal line, rho >= -d inside, rho == -d exactly on the medial axis,
//         and rho == 0 when the foot point is a sharp convex edge or vertex.
// ============================================================================

#ifndef SDF_GEMS_HAS_MD_MEDIAL_RADIUS
#define SDF_GEMS_HAS_MD_MEDIAL_RADIUS

// Value returned for an infinite medial radius (exterior of a convex region).
#define SDF_GEMS_MD_INFINITY 1e20

// ---------------------------------------------------------------------------
// Operators on a (d, rho) pair
// ---------------------------------------------------------------------------

// Offset by t (t > 0 grows the body, t < 0 erodes it). The medial axis does not move, so rho shifts with d.
// Growing is exact for convex bodies (in general for t below the exterior reach).
// Eroding by e = -t is exact where rho >= e or d <= -e (see mdOffsetIsExact); elsewhere the returned d is
// still a conservative (under-estimating) bound of the eroded body and the returned rho is clamped to 0.
vec2 mdOffset( vec2 m, float t )
{
    return vec2(m.x - t, max(m.y + t, 0.0));
}

// True iff d + e is the EXACT signed distance to the inner parallel body {d <= -e} at this point (e >= 0).
bool mdOffsetIsExact( vec2 m, float e )
{
    return m.x <= -e || m.y >= e;
}

// Local thickness at a surface point: diameter of the largest inscribed ball touching there.
// It is a lower bound of the chord length along the inward normal.
float mdThickness( vec2 m )
{
    return 2.0 * m.y;
}

// Distance from p to the medial axis measured ALONG ITS OWN NORMAL LINE (>= 0 inside, 0 on the medial axis).
// This is an UPPER bound of the Euclidean distance to the medial axis. For sphere tracing multiply by the
// measured per-primitive factor in the docs (1 for sphere/capsule/torus, 1/sqrt(2) for box and cylinder).
float mdSkeleton( vec2 m )
{
    return m.x + m.y;
}

// ---------------------------------------------------------------------------
// Primitives
// ---------------------------------------------------------------------------

// Sphere of radius r. Medial axis: the centre.
vec2 mdSphere( vec3 p, float r )
{
    return vec2(length(p) - r, r);
}

// Box with half-extents b. Medial axis: pieces of the bisector planes of the faces. 0 SFU for rho.
// With g = distances to the three face-plane pairs and i the nearest pair: rho = min(g_j, g_k, b_i).
vec2 mdBox( vec3 p, vec3 b )
{
    vec3  q = abs(p) - b;
    float m = max(q.x, max(q.y, q.z));
    float d = (max(min(q.x, q.y), min(max(q.x, q.y), q.z)) <= 0.0) ? m : length(max(q, 0.0));
    vec3  o = (q.x >= m) ? vec3(b.x, -q.y, -q.z)
            : (q.y >= m) ? vec3(-q.x, b.y, -q.z)
            :              vec3(-q.x, -q.y, b.z);
    return vec2(d, max(min(o.x, min(o.y, o.z)), 0.0));
}

// Round box: outer half-extents b, corner radius r (same convention as sdRoundBox), 0 <= r <= min(b).
vec2 mdRoundBox( vec3 p, vec3 b, float r )
{
    return mdOffset(mdBox(p, b - r), r);
}

// Capsule from a to b with radius r. Medial axis: the segment.
vec2 mdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
    vec3  pa = p - a, ba = b - a;
    float h  = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return vec2(length(pa - ba * h) - r, r);
}

// Capped cylinder along y, half-height h, radius r (same convention as sdCappedCylinder).
// Medial axis: two 45-degree cones from the rims plus an axis segment (h > r) or a central disk (h < r).
vec2 mdCylinder( vec3 p, float h, float r )
{
    float w = length(p.xz);
    float u = w - r, v = abs(p.y) - h;
    float d = (max(u, v) <= 0.0) ? max(u, v)
            : (v <= 0.0) ? u
            : (u <= 0.0) ? v
            : length(vec2(u, v));
    float rho = (u >= v) ? min(-v, r) : min(-u, h);
    return vec2(d, max(rho, 0.0));
}

// Ring torus in the xz-plane, t = (major R, minor r), R > r. Interior medial axis: the core circle.
vec2 mdTorus( vec3 p, vec2 t )
{
    return vec2(length(vec2(length(p.xz) - t.x, p.y)) - t.y, t.y);
}

// Ring torus, EXTERIOR medial radius: distance from the foot point along the OUTWARD normal to the exterior
// medial axis (the y axis). Minimum R - r (the exterior reach); SDF_GEMS_MD_INFINITY on the outer half.
vec2 mdTorusExterior( vec3 p, vec2 t )
{
    float w = length(p.xz);
    vec2  k = vec2(w - t.x, p.y);
    float l = length(k);
    float rho = (k.x < 0.0) ? min(max(t.x * l / (-k.x), t.x) - t.y, SDF_GEMS_MD_INFINITY) : SDF_GEMS_MD_INFINITY;
    return vec2(l - t.y, rho);
}

// Finite cone, apex at the origin, base at y = -h, c = (sin, cos) of the half-angle (same convention as sdCone).
// Medial axis: the axis from the apex to the in-centre, plus the base/mantle bisector cone from the rim.
vec2 mdCone( vec3 p, vec2 c, float h )
{
    float x2 = dot(p.xz, p.xz), w = sqrt(x2);
    float ds = c.y * w + c.x * p.y;    // signed distance to the mantle's supporting cone
    float db = -p.y - h;               // signed distance to the base plane
    float l  = c.x * w - c.y * p.y;    // slant coordinate of the foot on the mantle (0 at apex, h / c.y at rim)
    float rb = w - h * c.x / c.y;      // radial offset from the rim

    float d;
    if( max(ds, db) <= 0.0 )                        d = max(ds, db);
    else if( db >= 0.0 && rb <= 0.0 )               d = db;
    else if( l <= 0.0 )                             d = sqrt(x2 + p.y * p.y);
    else if( c.x * rb + c.y * db >= 0.0 )           d = length(vec2(rb, db));
    else                                            d = ds;

    float rho = (ds >= db) ? min(l * c.x / c.y, (h - l * c.y) / (1.0 + c.x))
                           : (h * c.x - w * c.y) / (1.0 + c.x);
    return vec2(d, max(rho, 0.0));
}

// Spheroid (ellipsoid of revolution about y): ab = (equatorial radius, polar radius), as sdSpheroid.
// Medial axis: the focal disk (oblate) or the focal segment (prolate). With m <= M the two semi-axes and
// (c, s) the direction cosines of the foot point (m c, M s):  rho = (m / M) sqrt(m^2 s^2 + M^2 c^2),
// ranging from m^2 / M (the inner reach, at the ends of the long axis) to m.
// The foot point comes from the same closed-form quartic solve as sdSpheroid.
vec2 mdSpheroid( vec3 p, vec2 ab )
{
    vec2 q = vec2( length(p.xz), abs(p.y) );
    if( ab.x > ab.y ) { q = q.yx; ab = ab.yx; }

    float a2 = ab.x * ab.x, b2 = ab.y * ab.y;
    float l  = b2 - a2, e = 0.5 * l;
    float x  = ab.x * q.x, y = ab.y * q.y;
    float x2 = x * x, y2 = y * y, r2 = x2 + y2;

    float C   = (r2 - l * l) / 3.0;
    float C2  = C * C, C3 = C2 * C;
    float lxy = l * x * y, l2x2y2 = lxy * lxy;
    float D   = C3 + l2x2y2, Q = C3 + 2.0 * l2x2y2;

    float v, vm;
    if( D < 0.0 )
    {
        float h   = atan( 2.0 * lxy * sqrt(-D), -Q ) / 3.0;
        float s   = cos(h), t = sin(h) * 1.732050808;
        float oms = (t * t) / (3.0 * (1.0 + s));
        float rz1 = sqrt( max( l * l + 2.0 * C * oms, 0.0 ) );
        float rz3 = sqrt( max( r2 - C * (oms + t), 0.0 ) );
        float rz2 = sqrt( max( r2 - C * (oms - t), 0.0 ) );
        v  = 0.5 * (rz1 + rz3 + ((x2 >= y2) ? rz2 : -rz2));
        vm = (x2 < y2) ? (-C) * (t / max(rz3 + rz2, 1e-30) - oms / max(rz1 + l, 1e-30)) : (v - e);
    }
    else
    {
        float H   = 2.0 * lxy * sqrt(D);
        float S   = pow( max(Q + H, 0.0), 1.0 / 3.0 );
        float U   = (S > 0.0) ? (C2 / S) : 0.0;
        float Zx  = r2 - C + 0.5 * (S + U);
        float dSU = S * S + C2 + U * U;
        float Zy  = 0.866025404 * ((dSU > 0.0) ? (2.0 * H / dSU) : 0.0);
        float Zm  = length( vec2(Zx, Zy) );
        v  = 0.5 * (sqrt( 2.0 * (Zm + Zx) ) + ((Zm > 0.0) ? l * (x2 - y2) / Zm : 0.0));
        vm = max( v - e, 0.0 );
    }

    vec2  w    = vec2( x * (v + e), y * max(vm, 0.0) );
    float wl   = length(w);
    float si_u = (l > 0.0) ? min( y / l, 1.0 ) : 0.0;
    vec2  cs   = (wl > 0.0) ? (w / wl) : vec2( sqrt(max(1.0 - si_u * si_u, 0.0)), si_u );
    float d    = length( q - ab * cs );
    d = (q.x * q.x * b2 + q.y * q.y * a2 >= a2 * b2) ? d : -d;

    float rho = (ab.x / ab.y) * length(vec2(ab.x * cs.y, ab.y * cs.x));
    return vec2(d, rho);
}

// Triaxial ellipsoid with semi-axes r: medial radius at a SURFACE point x (the foot point).
// The medial axis is the focal elliptical disk in the plane orthogonal to the shortest axis c, and
// rho(x) = c^2 * |(x / r^2)|.
float mdEllipsoidFoot( vec3 x, vec3 r )
{
    float c = min(r.x, min(r.y, r.z));
    return c * c * length(x / (r * r));
}

// Triaxial ellipsoid from a distance + unit gradient dg = (d, grad d) at p, e.g. sdgEllipsoidExact(p, r).
vec2 mdEllipsoidFromSdg( vec4 dg, vec3 p, vec3 r )
{
    return vec2(dg.x, mdEllipsoidFoot(p - dg.x * dg.yzw, r));
}

#endif // SDF_GEMS_HAS_MD_MEDIAL_RADIUS
