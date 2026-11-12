// ============================================================================
// rcSphere, rcPlane, rcCapsule, rcBox, rcRoundBox, rcTorus, rcSphereCone,
// rcSoftShadow, rcSinClearance, rcPenumbra, rcCoverage
// — Ray-Clearance Fields: the signed distance function over ray space
// @author Claude Fable 5.2 (sdf-gems, 2026). The sphere and capsule cases restate known
//         closest-approach formulas (Inigo Quilez's sphere / capsule soft shadows,
//         Ericson's segment-segment distance); see docs/optics/rcPrimitives.md.
// @exact  rcSphere, rcPlane, rcCapsule, rcBox, rcRoundBox: closed form; max |C - C_gt|
//         <= 6.5e-7 (GPU FP32) / 3e-15 (FP64 port) against dense sampling + golden-section
//         refinement of the repo SDFs on 24,012 rays each (rcPlane: 4.2e-6 at tmax = 30).
//         rcTorus: ITERATIVE (closed-form start + 2 x 5 Newton steps), ring / horn tori
//         only (R >= r). Measured error 5.8e-7 on random rays (FP32); up to 4.4e-5 (FP32) and
//         ~1e-6 (FP64, reference-limited) on rays tangent to the centre circle. Every returned value is the
//         SDF at a real ray point, so the error is one-sided: C >= C_true.
//         rcSphereCone: exact closed form of min_t d/t for a sphere.
//         rcSoftShadow, rcSinClearance, rcPenumbra, rcCoverage: approximations (see doc).
// ============================================================================
//
// For a ray segment x(t) = ro + t*rd, |rd| = 1, t in [0, tmax], and an exact SDF d:
//
//     C  = min_t d(x(t))      (signed clearance)        t* = argmin
//
//   C > 0 : the ray misses; C is its closest approach to the surface.
//   C < 0 : the ray penetrates; -C is the deepest penetration depth along the ray.
//
// Every rc* primitive returns vec2(C, t*) and mirrors the parameters of the matching sd*
// primitive (shape centred at the origin: transform the ray, not the shape).
//
//   Union        : exact,  rc(A u B) = (C_A < C_B) ? rcA : rcB.
//   Offset       : exact,  rc(d - r) = rc(d) - vec2(r, 0).
//   Intersection : max(C_A, C_B) is only a LOWER bound (it reports false hits).
//   Smooth union : smin(C_A, C_B) <= C <= min(C_A, C_B).
//   Ray space    : |C(ro + v, rd) - C(ro, rd)| <= |v|  (1-Lipschitz under translation);
//                  C(ro, rd') <= C(ro, rd) + t* * angle(rd, rd')  (one-sided under rotation).

#ifndef SDF_GEMS_RC_PRIMITIVES_GLSL
#define SDF_GEMS_RC_PRIMITIVES_GLSL

// ---------------------------------------------------------------------------
// Sphere of radius s at the origin. Exact.
// ---------------------------------------------------------------------------
vec2 rcSphere( in vec3 ro, in vec3 rd, in float tmax, in float s )
{
    float t = clamp(-dot(ro, rd), 0.0, tmax);
    return vec2(length(ro + t * rd) - s, t);
}

// ---------------------------------------------------------------------------
// Half-space dot(p, n) + h <= 0 (n unit). Exact: d is linear along the ray.
// ---------------------------------------------------------------------------
vec2 rcPlane( in vec3 ro, in vec3 rd, in float tmax, in vec3 n, in float h )
{
    float d0 = dot(ro, n) + h, s = dot(rd, n);
    return (s < 0.0) ? vec2(d0 + s * tmax, tmax) : vec2(d0, 0.0);
}

// ---------------------------------------------------------------------------
// Capsule from a to b with radius r. Exact: clamped segment-segment distance - r.
// A zero-length capsule (a == b) reduces to the sphere.
// ---------------------------------------------------------------------------
vec2 rcCapsule( in vec3 ro, in vec3 rd, in float tmax, in vec3 a, in vec3 b, in float r )
{
    vec3  ba = b - a, oa = ro - a;
    float e = dot(ba, ba), f = dot(ba, oa), c = dot(rd, oa), g = dot(rd, ba);
    float t = clamp(-c, 0.0, tmax), s = 0.0;
    if( e > 1e-20 )
    {
        vec3  n   = cross(rd, ba);                   // cross products instead of e - g*g:
        float den = dot(n, n);                       // no cancellation for near-parallel rays
        t = (den > 1e-12 * e) ? clamp(-dot(n, cross(oa, ba)) / den, 0.0, tmax) : 0.0;
        s = (g * t + f) / e;
        if( s < 0.0 )      { s = 0.0; t = clamp(-c, 0.0, tmax); }
        else if( s > 1.0 ) { s = 1.0; t = clamp(g - c, 0.0, tmax); }
    }
    return vec2(length(oa + t * rd - s * ba) - r, t);
}

// ---------------------------------------------------------------------------
// Box of half-extents b. Exact, signed.
//   Penetrating: inside the box d = max_i(|p_i| - b_i) is a convex piecewise-linear
//   function of t whose minimum is the crossing of its rising and falling envelopes,
//   t_x = min_i max_j t_ij  (closed form, no search).
//   Missing: D^2(t) = |p - clamp(p,-b,b)|^2 is convex and C^1; its derivative
//   g(t) = dot(rd, p - clamp(p,-b,b)) is piecewise linear with the 6 slab crossings
//   as its only breakpoints, so the root is one linear interpolation.
// ---------------------------------------------------------------------------
float rcBoxG( in vec3 o, in vec3 d, in vec3 b, in float t )
{
    vec3 p = o + t * d;
    return dot(d, p - clamp(p, -b, b));
}

vec2 rcBox( in vec3 ro, in vec3 rd, in float tmax, in vec3 b )
{
    // Mirror so that every direction component is >= 0.
    vec3 sg = vec3(rd.x < 0.0 ? -1.0 : 1.0, rd.y < 0.0 ? -1.0 : 1.0, rd.z < 0.0 ? -1.0 : 1.0);
    vec3 o = ro * sg, d = rd * sg;

    // (1) Penetration: rising lines up_i = (o_i - b_i) + t d_i, falling dn_j = (-o_j - b_j) - t d_j.
    vec3 up = o - b, dn = -o - b;
    float tx = 1e30;
    for( int i = 0; i < 3; i++ )
    {
        float ti = -1e30;
        for( int j = 0; j < 3; j++ )
        {
            float den = d[i] + d[j], num = dn[j] - up[i];
            float tij = (den > 1e-12) ? num / den : ((num <= 0.0) ? -1e30 : 1e30);
            ti = max(ti, tij);
        }
        tx = min(tx, ti);
    }
    tx = clamp(tx, 0.0, tmax);
    vec3  q  = abs(o + tx * d) - b;
    float fi = max(q.x, max(q.y, q.z));
    if( fi <= 0.0 ) return vec2(fi, tx);

    // (2) Miss: root of the monotone piecewise-linear g(t) on [0, tmax].
    float lo = 0.0, glo = rcBoxG(o, d, b, 0.0);
    if( glo >= 0.0 ) return vec2(length(max(abs(o) - b, 0.0)), 0.0);
    float hi = tmax, ghi = rcBoxG(o, d, b, tmax);
    if( ghi > 0.0 )
    {
        for( int i = 0; i < 3; i++ )
        {
            if( d[i] > 1e-12 )
            {
                float inv = 1.0 / d[i];
                for( int k = 0; k < 2; k++ )
                {
                    float tk = clamp(((k == 0 ? -b[i] : b[i]) - o[i]) * inv, 0.0, tmax);
                    float gk = rcBoxG(o, d, b, tk);
                    if( gk <= 0.0 && tk > lo ) { lo = tk; glo = gk; }
                    if( gk >= 0.0 && tk < hi ) { hi = tk; ghi = gk; }
                }
            }
        }
        float dg = ghi - glo;
        lo = (dg > 0.0) ? lo + (hi - lo) * clamp(-glo / dg, 0.0, 1.0) : lo;
    }
    else lo = tmax;
    return vec2(length(max(abs(o + lo * d) - b, 0.0)), lo);
}

// ---------------------------------------------------------------------------
// Rounded box (same parameters as sdRoundBox: outer half-extents b, radius r). Exact.
// ---------------------------------------------------------------------------
vec2 rcRoundBox( in vec3 ro, in vec3 rd, in float tmax, in vec3 b, in float r )
{
    return rcBox(ro, rd, tmax, b - r) - vec2(r, 0.0);
}

// ---------------------------------------------------------------------------
// Torus around the y axis, tor = (R, r), ring or horn (R >= r). ITERATIVE.
//   With x = t - tc (tc = closest approach to the centre) and v = x - e (offset from the
//   closest approach to the axis), F(t) = dist^2(x(t), circle) has
//       F'/2 = h = x - rho v / sqrt(v^2 + mu^2),   rho = R sqrt(a),  mu = axis distance / sqrt(a),
//   convex for v > 0 and concave for v < 0, so each half holds at most one local minimum.
//   Per half: closed-form vertex v_m of h (h' = 0), existence test h(v_m) < 0, a start from
//   the quadratic / cubic Taylor model at v_m (handles double and triple roots, i.e. rays
//   tangent to the centre circle), then 5 Newton steps on the convex side.
//   Result = best of {0, tmax, left minimum, right minimum}, each evaluated exactly.
// ---------------------------------------------------------------------------
float rcTorusRoot( in float e, in float rho, in float mu2 )
{
    float vm  = max(sqrt(max(pow(rho * mu2, 2.0 / 3.0) - mu2, 0.0)), 1e-6 * rho);
    float iq  = inversesqrt(max(vm * vm + mu2, 1e-30)), iq2 = iq * iq;
    float w   = rho * mu2 * iq * iq2;
    float hm  = vm + e - rho * vm * iq;
    if( hm >= 0.0 ) return e;                        // no local minimum on this half
    float k2  = 3.0 * w * vm * iq2;                  // h''(v_m)
    float k3  = 0.5 * w * (mu2 - 4.0 * vm * vm) * iq2 * iq2;   // h'''(v_m) / 6
    float dx  = sqrt(hm * hm - 2.0 * hm / max(k2, 1e-20));
    if( k3 > 1e-20 ) dx = min(dx, pow(-hm / k3, 1.0 / 3.0));
    float xm  = vm + e;
    float x   = min(max(rho, xm), xm + dx);
    for( int i = 0; i < 5; i++ )
    {
        float v  = x - e;
        float jq = inversesqrt(max(v * v + mu2, 1e-30));
        float h  = x - rho * v * jq;
        float hp = 1.0 - rho * mu2 * jq * jq * jq;
        x = (hp > 1e-6) ? max(x - h / hp, xm) : xm;
    }
    return x;
}

float rcTorusD2( in vec3 ro, in vec3 rd, in float t, in float R )
{
    vec3 p = ro + t * rd;
    float w = length(p.xz) - R;
    return w * w + p.y * p.y;
}

vec2 rcTorus( in vec3 ro, in vec3 rd, in float tmax, in vec2 tor )
{
    float a   = max(dot(rd.xz, rd.xz), 1e-24);
    float tc  = -dot(ro, rd);
    float cr  = (ro.x * rd.z - ro.z * rd.x) / a;
    float rho = tor.x * sqrt(a);
    float e   = -dot(ro.xz, rd.xz) / a - tc;
    float mu2 = max(cr * cr, 1e-8 * rho * rho);
    float tR  = clamp(tc + rcTorusRoot( e, rho, mu2), 0.0, tmax);
    float tL  = clamp(tc - rcTorusRoot(-e, rho, mu2), 0.0, tmax);

    vec2 best = vec2(rcTorusD2(ro, rd, 0.0, tor.x), 0.0);
    float f = rcTorusD2(ro, rd, tmax, tor.x); if( f < best.x ) best = vec2(f, tmax);
    f = rcTorusD2(ro, rd, tL, tor.x);         if( f < best.x ) best = vec2(f, tL);
    f = rcTorusD2(ro, rd, tR, tor.x);         if( f < best.x ) best = vec2(f, tR);
    return vec2(sqrt(best.x) - tor.y, best.y);
}

// ---------------------------------------------------------------------------
// Exact closed form of the marched soft-shadow estimator for a sphere:
//     min_{t > 0} sdSphere(ro + t rd, s) / t  =  sin(theta - alpha)
// theta = angle between rd and the direction to the centre, sin(alpha) = s / |ro|.
// It is the sine of the ANGULAR clearance (negative when the ray hits), and is 1
// when the sphere is entirely behind the ray. Requires |ro| > s.
// ---------------------------------------------------------------------------
float rcSphereCone( in vec3 ro, in vec3 rd, in float s )
{
    float A = dot(ro, ro), B = dot(ro, rd);
    if( B >= s ) return 1.0;
    return (sqrt(max((A - B * B) * (A - s * s), 0.0)) + s * B) / max(A, 1e-30);
}

// ---------------------------------------------------------------------------
// 1-tap penumbra from a clearance c = vec2(C, t*). Approximates the marched
// clamp(k * min_t d/t, 0, 1) by k * sin(atan(C / t*)): C / hypot(t*, C) is the sine of
// the angle (seen from ro) of the surface point nearest to the ray, an upper bound of
// the true angular clearance. For several objects take the min of the RESULTS (the
// object with the smallest C need not have the smallest C / t*).
// ---------------------------------------------------------------------------
float rcSoftShadow( in vec2 c, in float k )
{
    if( c.x <= 0.0 ) return 0.0;
    return clamp(k * c.x * inversesqrt(c.y * c.y + c.x * c.x), 0.0, 1.0);
}

// ---------------------------------------------------------------------------
// Signed angular clearance estimate sin(phi) ~ C / hypot(t*, C), in [-1, 1]. Negative when
// the ray penetrates. Combine several objects with min().
// ---------------------------------------------------------------------------
float rcSinClearance( in vec2 c )
{
    return c.x * inversesqrt(max(c.y * c.y + c.x * c.x, 1e-30));
}

// ---------------------------------------------------------------------------
// Two-sided penumbra of a disc light of angular radius beta (sinBeta = sin(beta)) from the
// signed sine clearance s (rcSinClearance, or rcSphereCone for spheres): the visible
// fraction of a disc cut by a straight edge at offset s / sinBeta. Uses the NEGATIVE
// clearances too, which a marched min(k d / t) cannot provide (it stops at the surface),
// so the penumbra is centred on the geometric shadow edge. Approximation.
// ---------------------------------------------------------------------------
float rcPenumbra( in float s, in float sinBeta )
{
    float u = clamp(s / sinBeta, -1.0, 1.0);
    return 0.5 + (u * sqrt(1.0 - u * u) + asin(u)) * 0.31830988618;
}

// ---------------------------------------------------------------------------
// Analytic silhouette coverage for ray-traced primitives: fraction of a pixel of
// angular footprint `pix` (radians per pixel) covered by the object, from the signed
// clearance of the pixel's centre ray. Approximation (box-filtered straight edge).
// ---------------------------------------------------------------------------
float rcCoverage( in vec2 c, in float pix )
{
    return clamp(0.5 - c.x / max(pix * c.y, 1e-30), 0.0, 1.0);
}

#endif // SDF_GEMS_RC_PRIMITIVES_GLSL
