// ============================================================================
// opLaguerreBoost — Lorentz boost of the eikonal graph (a Laguerre transformation of an SDF)
//   opLaguerreBoostInit / Point / Step / Finalize, sdgLaguerreBoostFinalize, opLaguerreBoostProbe,
//   sdLaguerreSphere, sdLaguerrePlane, sdLaguerreCapsule, sdLaguerreBox, sdLaguerreTorus (+ sdg variants)
// @author Claude Fable 5.2 (sdf-gems, 2026)
// @exact  The operator output d' is ALWAYS 1-Lipschitz with |grad d'| = 1 a.e., but it is the exact SDF of its
//         own zero set only where the null generator through the evaluated graph point reaches t' = 0 before
//         ending on the medial axis of the input shape (see docs/deform/opLaguerreBoost.md, "Findings").
//         sdLaguerreSphere / sdLaguerrePlane: Exact.
//         sdLaguerreCapsule: Exact (a round cone) iff r >= |beta| h; otherwise Conservative bound.
//         sdLaguerreTorus:   Exact (a ring Dupin cyclide) iff r >= |beta| R (up to the 8-tap Newton solve, |beta| <= 0.95);
//                            otherwise Conservative bound.
//         sdLaguerreBox:     Closed form of the operator output. Conservative bound in general (exact inside and
//                            in the regions listed in the doc); its zero set is NOT a polyhedron unless beta = 0.
//
// Convention: the graph {(x, t) : t = d(x)} of an SDF is a null hypersurface of R^{3,1} (signature +++-).
// Boost along the unit axis e with velocity beta in (-1, 1), gamma = 1 / sqrt(1 - beta^2):
//     x'_par = gamma (x_par - beta t),   x'_perp = x_perp,   t' = gamma (t - beta x_par).
// To evaluate d'(p): find the scalar s (= x_par of the pre-image) with
//     s = u + beta * d(perp + s e),   u = dot(p, e) / gamma,   perp = p - dot(p, e) e,
// a contraction with factor |beta| because d is 1-Lipschitz. Then d' = gamma (d - beta s) = d / gamma - gamma beta u.
// ============================================================================

#ifndef SDF_GEMS_HAS_OP_LAGUERRE_BOOST
#define SDF_GEMS_HAS_OP_LAGUERRE_BOOST

// ---------------------------------------------------------------------------
// 1. Generic operator, expressed as helpers around the taps of the inner SDF.
//
//      vec3  perp;
//      float u = opLaguerreBoostInit(p, e, beta, perp), s = u;
//      for (int i = 0; i < 4; i++) {                                   // Newton taps: see the doc table
//          vec4 dg = sdgInner(opLaguerreBoostPoint(s, perp, e));
//          s = opLaguerreBoostStep(s, u, beta, e, dg);
//      }
//      vec4 r = sdgLaguerreBoostFinalize(s, u, beta, e, sdgInner(opLaguerreBoostPoint(s, perp, e)));
// ---------------------------------------------------------------------------

// Returns u = dot(p, e) / gamma (also the initial guess for s) and the component of p orthogonal to e.
float opLaguerreBoostInit( in vec3 p, in vec3 e, in float beta, out vec3 perp )
{
    float x = dot(p, e);
    perp = p - x * e;
    return x * sqrt((1.0 - beta) * (1.0 + beta));
}

// Pre-image sample position for the current s.
vec3 opLaguerreBoostPoint( in float s, in vec3 perp, in vec3 e )
{
    return perp + s * e;
}

// Fixed-point step (needs only the distance): error shrinks by |beta| per tap.
float opLaguerreBoostStep( in float u, in float beta, in float d )
{
    return u + beta * d;
}

// Newton step (needs the gradient): F(s) = s - u - beta d(s), F' = 1 - beta d_s in [1 - |beta|, 1 + |beta|].
float opLaguerreBoostStep( in float s, in float u, in float beta, in vec3 e, in vec4 dg )
{
    return s - (s - u - beta * dg.x) / (1.0 - beta * dot(dg.yzw, e));
}

// Safeguarded Newton step for NON-CONVEX inner SDFs with |beta| > 1/3, where plain Newton can cycle or diverge
// (measured on a torus at |beta| >= 0.8). Every tap proves that the root lies between s - F / (1 - |beta|) and
// s - F / (1 + |beta|); `br` (initialise to vec2(-1e30, 1e30)) accumulates that bracket, and a Newton step that
// leaves it is replaced by the bracket midpoint.
float opLaguerreBoostStepSafe( in float s, in float u, in float beta, in vec3 e, in vec4 dg, inout vec2 br )
{
    float F  = s - u - beta * dg.x;
    float a  = s - F / (1.0 - abs(beta)), b = s - F / (1.0 + abs(beta));
    br = vec2(max(br.x, min(a, b)), min(br.y, max(a, b)));
    float sn = s - F / (1.0 - beta * dot(dg.yzw, e));
    return (sn >= br.x && sn <= br.y) ? sn : 0.5 * (br.x + br.y);
}

// Value from the last tap d = d(perp + s e). Written as d / gamma - gamma beta u (equal to gamma (d - beta s) at
// the root); with an unconverged s its error is bounded by the tap residual times the Doppler factor:
// |error| <= |s - u - beta d| * sqrt((1 + |beta|) / (1 - |beta|)).
float opLaguerreBoostFinalize( in float u, in float beta, in float d )
{
    float ig = sqrt((1.0 - beta) * (1.0 + beta));
    return d * ig - beta * u / ig;
}

// Value + gradient from the last tap dg = sdg(perp + s e). The value includes one free Newton correction
// (second-order accurate in the remaining residual). The gradient is the relativistic aberration of the normal:
//   n'_par = (n_par - beta) / (1 - beta n_par),   n'_perp = n_perp / (gamma (1 - beta n_par)),   |n'| = |n| = 1.
vec4 sdgLaguerreBoostFinalize( in float s, in float u, in float beta, in vec3 e, in vec4 dg )
{
    float ig = sqrt((1.0 - beta) * (1.0 + beta));
    float np = dot(dg.yzw, e);
    float k  = 1.0 / (1.0 - beta * np);
    float d  = dg.x - np * (s - u - beta * dg.x) * k;
    vec3  g  = (ig * k) * (dg.yzw - np * e) + ((np - beta) * k) * e;
    return vec4(d * ig - beta * u / ig, g);
}

// Exactness probe (one extra tap). x = perp + s e is the solved pre-image point and dg = sdg(x).
// Returns vec4(z, s0): d' is the exact distance at p iff the inner SDF satisfies d(z) == s0, i.e. the normal line
// through x is still a generator of the input SDF where it crosses the new zero level t' = 0.
vec4 opLaguerreBoostProbe( in vec3 x, in vec3 e, in float beta, in vec4 dg )
{
    vec3  y  = x - dg.x * dg.yzw;
    float s0 = beta * dot(y, e) / (1.0 - beta * dot(dg.yzw, e));
    return vec4(y + s0 * dg.yzw, s0);
}

// ---------------------------------------------------------------------------
// 2. Closed forms: sphere -> sphere, plane -> plane (exact).
// ---------------------------------------------------------------------------

// Sphere (centre c, radius r): the cone apex (c, -r) is boosted to centre c + ((gamma - 1) c_par + gamma beta r) e
// and radius gamma (r + beta c_par). A negative radius is a legitimate oriented sphere: d' > 0 everywhere.
float sdLaguerreSphere( in vec3 p, in vec3 c, in float r, in vec3 e, in float beta )
{
    float ga = inversesqrt((1.0 - beta) * (1.0 + beta));
    float cp = dot(c, e);
    return length(p - c - (ga * (cp + beta * r) - cp) * e) - ga * (r + beta * cp);
}

// Plane dot(n, x) - h (|n| = 1): normal is aberrated, offset becomes h / (gamma (1 - beta n_par)).
vec4 sdgLaguerrePlane( in vec3 p, in vec3 n, in float h, in vec3 e, in float beta )
{
    float ig = sqrt((1.0 - beta) * (1.0 + beta));
    float np = dot(n, e);
    float k  = 1.0 / (1.0 - beta * np);
    vec3  g  = (ig * k) * (n - np * e) + ((np - beta) * k) * e;
    return vec4(dot(p, g) - h * ig * k, g);
}

float sdLaguerrePlane( in vec3 p, in vec3 n, in float h, in vec3 e, in float beta )
{
    return sdgLaguerrePlane(p, n, h, e, beta).x;
}

// ---------------------------------------------------------------------------
// 3. Capsule on the boost axis -> round cone.
//    Capsule: segment (-h, 0, 0)..(h, 0, 0), radius r, boosted along +x. The medial spheres (s, r) map to
//    centres gamma (s + beta r) with radii gamma (r + beta s): centre and radius are both linear in s with
//    d(radius) / d(centre) = beta, i.e. a round cone whose half-angle is asin(beta).
//    Exact iff both radii are >= 0, i.e. r >= |beta| h. Otherwise the zero set is a cone with a sharp apex
//    plus a spherical cap, and the value is a conservative bound outside near the apex.
// ---------------------------------------------------------------------------
vec4 sdgLaguerreCapsule( in vec3 p, in float h, in float r, in float beta )
{
    float ga = inversesqrt((1.0 - beta) * (1.0 + beta));
    float c1 = ga * (beta * r - h), r1 = ga * (r - beta * h);   // image of the sphere at s = -h
    float c2 = ga * (beta * r + h), r2 = ga * (r + beta * h);   // image of the sphere at s = +h
    float x  = p.x - c1, z = p.x - c2;
    float l2 = dot(p.yz, p.yz);
    float a2 = (1.0 - beta) * (1.0 + beta);
    float k  = beta * abs(beta) * l2;                           // radius grows along +x at rate beta
    vec2  rd = (l2 > 1e-30) ? p.yz * inversesqrt(l2) : vec2(1.0, 0.0);
    if (a2 * x * abs(x) < -k)
    {
        float l = sqrt(l2 + x * x);
        return vec4(l - r1, (l > 1e-20) ? vec3(x, p.yz) / l : vec3(-1.0, 0.0, 0.0));
    }
    if (a2 * z * abs(z) > -k)
    {
        float l = sqrt(l2 + z * z);
        return vec4(l - r2, (l > 1e-20) ? vec3(z, p.yz) / l : vec3(1.0, 0.0, 0.0));
    }
    float ca = sqrt(a2);
    return vec4(ca * sqrt(l2) - beta * x - r1, -beta, ca * rd);
}

float sdLaguerreCapsule( in vec3 p, in float h, in float r, in float beta )
{
    return sdgLaguerreCapsule(p, h, r, beta).x;
}

// ---------------------------------------------------------------------------
// 4. Box boosted along +x: closed form of the operator output (0 iterations, 1-2 sqrt + 1 inversesqrt).
//    Along the line (s, y, z) the box SDF is constant, then linear or hyperbolic in s, so the fixed-point
//    equation is linear or quadratic. Always 1-Lipschitz with the correct zero set and sign; see the doc for
//    where it is the exact distance (faces with positive image offset, rounded +beta side) and where it is
//    only a bound (outside the sharp edges / vertices on the -beta side).
// ---------------------------------------------------------------------------
vec4 sdgLaguerreBox( in vec3 p, in vec3 b, in float beta )
{
    float a2 = (1.0 - beta) * (1.0 + beta);
    float ga = inversesqrt(a2);
    float u  = p.x * a2 * ga;
    vec2  q  = abs(p.yz) - b.yz;
    vec2  qp = max(q, 0.0);
    float m  = max(q.x, q.y);
    float e2 = dot(qp, qp);
    float cm = (m < 0.0) ? max(m, -b.x) : sqrt(e2);            // value of d on the central plateau
    float w0 = b.x + min(cm, 0.0);                              // plateau is |s| <= w0
    float t  = u + beta * cm;
    float s, d;
    vec3  n;                                                    // gradient of the box SDF at the pre-image
    if (abs(t) <= w0)
    {
        s = t; d = cm;
        n = (m < 0.0) ? ((q.x > q.y) ? vec3(0.0, 1.0, 0.0) : vec3(0.0, 0.0, 1.0))
                      : vec3(0.0, qp * inversesqrt(max(e2, 1e-30)));
    }
    else
    {
        float sg = (t >= 0.0) ? 1.0 : -1.0;
        float w  = sg * u - b.x;
        if (m < 0.0)
        {
            d = w / (1.0 - sg * beta);                          // d = |s| - b.x
            n = vec3(1.0, 0.0, 0.0);
        }
        else
        {
            float a = (w + sg * beta * sqrt(w * w + a2 * e2)) * ga * ga;   // a = |s| - b.x > 0
            d = sqrt(a * a + e2);
            n = vec3(a, qp) / d;
        }
        s = u + beta * d;
        n.x *= sg;
    }
    n.yz *= vec2((p.y >= 0.0) ? 1.0 : -1.0, (p.z >= 0.0) ? 1.0 : -1.0);
    float k = 1.0 / (1.0 - beta * n.x);
    return vec4(ga * (d - beta * s), (n.x - beta) * k, n.yz * (k * a2 * ga));
}

float sdLaguerreBox( in vec3 p, in vec3 b, in float beta )
{
    return sdgLaguerreBox(p, b, beta).x;
}

// ---------------------------------------------------------------------------
// 5. Torus boosted along a diameter -> Dupin cyclide (iterative instance of the generic operator).
//    Torus: axis z, major radius R, minor radius r, boosted along +x. The medial circle spheres map to spheres
//    centred on the ellipse (gamma (R cos t + beta r), R sin t, 0) with radii gamma (r + beta R cos t): a Dupin
//    cyclide whose focal ellipse has eccentricity exactly |beta|. Ring cyclide (exact) iff r >= |beta| R.
//    7 safeguarded Newton taps + 1 final tap (verified for |beta| <= 0.95); measured error in the doc.
// ---------------------------------------------------------------------------
vec4 sdgLaguerreTorusInner( in vec3 x, in float R, in float r )
{
    float l  = length(x.xy);
    vec2  w  = vec2(l - R, x.z);
    float lw = length(w);
    w = (lw > 1e-20) ? w / lw : vec2(1.0, 0.0);
    vec2  rd = (l > 1e-20) ? x.xy / l : vec2(1.0, 0.0);
    return vec4(lw - r, w.x * rd, w.y);
}

vec4 sdgLaguerreTorus( in vec3 p, in float R, in float r, in float beta )
{
    const vec3 e = vec3(1.0, 0.0, 0.0);
    vec3  perp;
    float u = opLaguerreBoostInit(p, e, beta, perp), s = u;
    vec2  br = vec2(-1e30, 1e30);
    for (int i = 0; i < 7; i++)
    {
        vec4 dg = sdgLaguerreTorusInner(opLaguerreBoostPoint(s, perp, e), R, r);
        s = opLaguerreBoostStepSafe(s, u, beta, e, dg, br);
    }
    return sdgLaguerreBoostFinalize(s, u, beta, e, sdgLaguerreTorusInner(opLaguerreBoostPoint(s, perp, e), R, r));
}

float sdLaguerreTorus( in vec3 p, in float R, in float r, in float beta )
{
    return sdgLaguerreTorus(p, R, r, beta).x;
}

#endif // SDF_GEMS_HAS_OP_LAGUERRE_BOOST
