// ============================================================================
// keRopelength — Gonzalez-Maddocks Ropelength Thickness, Tangent-Sphere Radii,
//                Mobius-Invariant Pair Potentials & Self-Contact Link SDFs
//   kePointTangentRadius, keSurfaceTangentRadius, keMobiusPairEnergy,
//   keHopfLinkInjectivityRadius, keTwoCircleLinkThickness,
//   sdTightHopfLink3D, sdgTightHopfLink3D, keApolloniusContactTube, sdSwollenHopfLink3D
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  sdTightHopfLink3D / sdgTightHopfLink3D: Exact Euclidean SDF and analytical unit gradient
//         everywhere in R^3 (both interior and exterior) for |eccentricity| < 1, because clamping
//         the tube radius to the ropelength injectivity bound r_max prevents interior overlap.
//         keApolloniusContactTube / sdSwollenHopfLink3D: Globally 1-Lipschitz conservative bound
//         (proven via 1 / (1 + |alpha|) normalization; scaling by 1 / sqrt(1 + alpha^2) fails
//         1-Lipschitz where gradients oppose).
// ============================================================================

#ifndef SDF_GEMS_HAS_KE_ROPELENGTH
#define SDF_GEMS_HAS_KE_ROPELENGTH

// ---------------------------------------------------------------------------
// 1. Gonzalez-Maddocks Point-Tangent Circle Radius & Surface-Tangent Sphere Radius
// ---------------------------------------------------------------------------

// Unique circle radius R_pt(x, t; y) through x and y tangent to unit curve tangent t at x
// (Gonzalez & Maddocks 1999). Uses length(v - dot(v, t) * t) rather than sqrt(|v|^2 - (v.t)^2)
// to avoid FP32 catastrophic cancellation when y - x is nearly parallel to t.
float kePointTangentRadius( in vec3 x, in vec3 t, in vec3 y )
{
    vec3  v  = y - x;
    float v2 = dot(v, v);
    vec3  vp = v - dot(v, t) * t;
    float l  = length(vp);
    return (l > 1e-20) ? (0.5 * v2 / l) : 1e20;
}

// Signed radius R_tan(p; x, n) of the sphere tangent to the oriented surface element (x, n)
// passing through p: R_tan = |p - x|^2 / (2 dot(p - x, n)).
// Positive when p lies in the +n half-space, negative in the -n half-space.
// The medial radius / reach at x is inf_{y in Sigma \ {x}} |R_tan(y; x, n)|.
float keSurfaceTangentRadius( in vec3 p, in vec3 x, in vec3 n )
{
    vec3  v   = p - x;
    float v2  = dot(v, v);
    float vn2 = 2.0 * dot(v, n);
    return (abs(vn2) > 1e-20) ? (v2 / vn2) : ((vn2 >= 0.0) ? 1e20 : -1e20);
}

// ---------------------------------------------------------------------------
// 2. Mobius-Invariant Pair Potential (Kusner-Sullivan / Freedman-He-Wang Conformal Kernel)
// ---------------------------------------------------------------------------

// Given two oriented curve elements (x, tx) and (y, ty) with |tx| = |ty| = 1:
// The circle C(x, tx, y) through x and y tangent to tx at x has tangent at y equal to the
// chord reflection t_tilde = 2 (tx . w) w - tx, where w = (y - x) / |y - x|.
// Returns vec2(K_mobius, one_minus_cos_phi):
//   .x = (1 - cos phi) / |y - x|^2  (conformal 2-form density; vanishes identically on circles)
//   .y = 1 - cos phi                (strictly invariant under the 3D Mobius group PO(4,1))
vec2 keMobiusPairEnergy( in vec3 x, in vec3 tx, in vec3 y, in vec3 ty )
{
    vec3  v   = y - x;
    float r2  = dot(v, v);
    if (r2 <= 1e-24) return vec2(0.0);
    float ir2 = 1.0 / r2;
    float vx  = dot(v, tx);
    float vy  = dot(v, ty);
    float omc = max(1.0 + dot(tx, ty) - 2.0 * vx * vy * ir2, 0.0);
    return vec2(omc * ir2, omc);
}

// ---------------------------------------------------------------------------
// 3. Closed-Form Two-Circle Link Thickness & Ropelength Injectivity Radius
// ---------------------------------------------------------------------------

// Exact maximum embedded tube radius r_max for two orthogonal circles:
//   C1 in the xy-plane, centered at (0, 0, 0), radius R1
//   C2 in the xz-plane, centered at (offset, 0, 0), radius R2
// Because dist(q, C1)^2 = R1^2 + R2^2 - a^2 + 2 a x - 2 R1 |x| is concave piecewise-linear
// in x in [a - R2, a + R2], the minimum mutual distance is ALWAYS attained at the x-axis
// pierce points x = a - R2 and x = a + R2.
// Combining mutual strand contact (0.5 * d_min) with local focal curvature bounds (R1, R2):
float keHopfLinkInjectivityRadius( in float R1, in float R2, in float offset )
{
    float a   = abs(offset);
    float dm1 = abs(abs(a - R2) - R1);
    float dm2 = abs(a + R2 - R1);
    float rMutual = 0.5 * min(dm1, dm2);
    return min(min(R1, R2), rMutual);
}

// Returns vec3(r_max, r_mutual, r_focal) for two orthogonal circles (R1, R2, offset):
//   .x = exact global radius of curvature (ropelength injectivity radius r_max)
//   .y = mutual strand contact bound 0.5 * dist(C1, C2)
//   .z = focal self-contact bound min(R1, R2)
vec3 keTwoCircleLinkThickness( in float R1, in float R2, in float offset )
{
    float a       = abs(offset);
    float rMutual = 0.5 * min(abs(abs(a - R2) - R1), abs(a + R2 - R1));
    float rFocal  = min(R1, R2);
    return vec3(min(rMutual, rFocal), rMutual, rFocal);
}

// ---------------------------------------------------------------------------
// 4. Exact Tight Hopf Link SDF & Analytical Gradient (sdTightHopfLink3D, sdgTightHopfLink3D)
// ---------------------------------------------------------------------------

// Helper: exact SDF + unit gradient of a circle tube in xy-plane (radius R, tube r) at q
vec4 keCircleTubeXY( in vec3 q, in float R, in float r )
{
    float lxy = length(q.xy);
    vec2  rad = (lxy > 1e-20) ? (q.xy / lxy) : vec2(1.0, 0.0);
    vec2  w   = vec2(lxy - R, q.z);
    float lw  = length(w);
    vec2  wn  = (lw > 1e-20) ? (w / lw) : vec2(1.0, 0.0);
    return vec4(lw - r, wn.x * rad.x, wn.x * rad.y, wn.y);
}

// Exact Euclidean SDF + unit gradient of a two-circle Hopf link centered at the origin:
//   R            : circle major radius (R > 0)
//   eccentricity : dimensionless offset parameter e in (-1, 1); center separation is a = (1 + e) R.
//                  At e = 0 (the ideal Gehring link), a = R and r_max = 0.5 * R (ropelength 8*pi).
//                  For general e in (-1, 1), r_max = 0.5 * R * (1 - |e|).
// Because r = r_max guarantees disjoint interiors, min(d1, d2) is the EXACT Euclidean SDF everywhere
// in R^3 (both interior and exterior).
vec4 sdgTightHopfLink3D( in vec3 p, in float R, in float eccentricity )
{
    float e    = clamp(eccentricity, -0.999, 0.999);
    float a    = (1.0 + e) * R;
    float rMax = 0.5 * R * (1.0 - abs(e));
    vec3  p1   = p - vec3(-0.5 * a, 0.0, 0.0);
    vec3  p2   = p - vec3( 0.5 * a, 0.0, 0.0);
    vec4  dg1  = keCircleTubeXY(p1, R, rMax);
    vec4  dg2  = keCircleTubeXY(vec3(p2.x, p2.z, p2.y), R, rMax);
    dg2.yzw    = vec3(dg2.y, dg2.w, dg2.z);
    return (dg1.x <= dg2.x) ? dg1 : dg2;
}

float sdTightHopfLink3D( in vec3 p, in float R, in float eccentricity )
{
    return sdgTightHopfLink3D(p, R, eccentricity).x;
}

// ---------------------------------------------------------------------------
// 5. 1-Lipschitz Apollonius / Bisector Contact-Limited Variable-Radius Tube
// ---------------------------------------------------------------------------

// Given two 1-Lipschitz core distances d1(p), d2(p) >= 0 (or signed SDFs), a maximum target radius r0,
// and contact saturation factor alpha in [0, 1]:
// Each tube swells up to r(p) = min(r0, alpha * 0.5 * (d1 + d2)) -> bounded by (d_near - alpha * d_far) / (1 + alpha).
// Dividing by (1 + |alpha|) (NOT sqrt(1 + alpha^2)) makes the contact-limited field globally 1-Lipschitz.
float keApolloniusContactTube( in float d1, in float d2, in float r0, in float alpha )
{
    float a  = clamp(alpha, 0.0, 1.0);
    float m  = min(d1, d2);
    float M  = max(d1, d2);
    float dc = (m - a * M) / (1.0 + a);
    return max(m - r0, dc);
}

// Self-avoiding variable-radius Hopf link that swells beyond r_max (for r0 > r_max) while clamping
// mutual contact along the exact Voronoi bisector d1 == d2 without interior overlap.
float sdSwollenHopfLink3D( in vec3 p, in float R, in float eccentricity, in float r0, in float alpha )
{
    float a  = (1.0 + clamp(eccentricity, -0.999, 0.999)) * R;
    vec3  p1 = p - vec3(-0.5 * a, 0.0, 0.0);
    vec3  p2 = p - vec3( 0.5 * a, 0.0, 0.0);
    float d1 = length(vec2(length(p1.xy) - R, p1.z));
    float d2 = length(vec2(length(p2.xz) - R, p2.y));
    return keApolloniusContactTube(d1, d2, r0, alpha);
}

#endif // SDF_GEMS_HAS_KE_ROPELENGTH
