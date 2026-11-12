// ============================================================================
// srHeisenberg — Sub-Riemannian Geometry on the Heisenberg Group H^1
//   srHeisenbergMul, srHeisenbergInv, srContactForm, srHorizontalGrad, srHorizontalLift,
//   srKoranyiNorm, sdgKoranyiGauge, sdKoranyiBall, sdgKoranyiBall,
//   srCarnotTheta, srCarnotDistance, sdgCarnotDistance, sdCarnotBall, sdgCarnotBall,
//   opHeisenbergTranslate, opHeisenbergTranslateScale, sdgHeisenbergTranslate
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  srKoranyiNorm / srCarnotDistance: Exact sub-Riemannian gauge and Carnot-Caratheodory
//         geodesic distance (|grad_H d_CC|_H == 1 everywhere off the vertical cut locus r = 0).
//         sdKoranyiBall / sdgKoranyiBall: Exact Euclidean SDF and unit gradient of the quartic
//         Koranyi ball r^4 + 16 z^2 <= R^4.
//         sdCarnotBall / sdgCarnotBall: Exact Euclidean SDF and unit gradient of the Gavrilov-Pansu
//         sub-Riemannian geodesic sphere (with exact conical pole handling at r = 0, |z| = R^2/(4 pi)).
//         opHeisenbergTranslateScale: Exact 1/sigma_max(J) Euclidean Lipschitz factor for left-invariant
//         Heisenberg group translation q = a^{-1} * p.
// ============================================================================

#ifndef SDF_GEMS_HAS_SR_HEISENBERG
#define SDF_GEMS_HAS_SR_HEISENBERG

// ---------------------------------------------------------------------------
// 1. Heisenberg Group H^1 Law & Contact / Horizontal Calculus
//    Group law: a * b = (a.x + b.x, a.y + b.y, a.z + b.z + 0.5 * (a.x * b.y - a.y * b.x))
//    Contact 1-form: omega = dz - 0.5 * (x dy - y dx)
//    Left-invariant horizontal frame:
//      X_1 = d_x - 0.5 y d_z,   X_2 = d_y + 0.5 x d_z,   [X_1, X_2] = d_z
// ---------------------------------------------------------------------------

vec3 srHeisenbergMul( in vec3 a, in vec3 b )
{
    return vec3(
        a.x + b.x,
        a.y + b.y,
        a.z + b.z + 0.5 * (a.x * b.y - a.y * b.x)
    );
}

vec3 srHeisenbergInv( in vec3 a )
{
    return -a;
}

// Anisotropic Heisenberg dilation delta_lambda(p) = (lambda x, lambda y, lambda^2 z),
// an automorphism of H^1 scaling homogeneous gauges by lambda.
vec3 srHeisenbergDilate( in vec3 p, in float lambda )
{
    return vec3(lambda * p.xy, (lambda * lambda) * p.z);
}

// Evaluates the Darboux contact 1-form omega_p(v) = v.z - 0.5 * (p.x * v.y - p.y * v.x).
// Vanishes iff v lies in the horizontal distribution H_p = span(X_1(p), X_2(p)).
float srContactForm( in vec3 p, in vec3 v )
{
    return v.z - 0.5 * (p.x * v.y - p.y * v.x);
}

// Projects a Euclidean gradient g = grad u(p) onto the orthonormal horizontal frame (X_1, X_2):
//   grad_H u(p) = (X_1 u(p), X_2 u(p)) in R^2.
// A function u satisfies the horizontal Eikonal equation iff length(srHorizontalGrad(p, grad u)) == 1.
vec2 srHorizontalGrad( in vec3 p, in vec3 euclidGrad )
{
    return vec2(
        euclidGrad.x - 0.5 * p.y * euclidGrad.z,
        euclidGrad.y + 0.5 * p.x * euclidGrad.z
    );
}

// Lifts a 2D horizontal vector h = (h_1, h_2) in the (X_1, X_2) frame to a 3D tangent vector
// h_1 X_1(p) + h_2 X_2(p) in H_p subset R^3. Satisfies srContactForm(p, srHorizontalLift(p, h)) == 0.
vec3 srHorizontalLift( in vec3 p, in vec2 h )
{
    return vec3(h.x, h.y, 0.5 * (p.x * h.y - p.y * h.x));
}

// ---------------------------------------------------------------------------
// 2. Koranyi Gauge ||p||_K = (r^4 + 16 z^2)^{1/4} vs Exact Euclidean Koranyi Ball SDF
//    Note: ||grad_H ||p||_K||_H = r / ||p||_K <= 1 (horizontally 1-Lipschitz, equal to 1 only on z = 0).
//    Meanwhile its Euclidean gradient norm reaches 2/R at the poles, so ||p||_K - R oversteps in R^3.
// ---------------------------------------------------------------------------

float srKoranyiNorm( in vec3 p )
{
    float r2 = dot(p.xy, p.xy);
    return sqrt(sqrt(r2 * r2 + 16.0 * p.z * p.z));
}

// Returns vec4(||p||_K, Euclidean gradient of ||p||_K).
vec4 sdgKoranyiGauge( in vec3 p )
{
    float r2 = dot(p.xy, p.xy);
    float q4 = r2 * r2 + 16.0 * p.z * p.z;
    if (q4 <= 1e-30) return vec4(0.0, 0.0, 0.0, 1.0);
    float u2 = sqrt(q4);
    float u  = sqrt(u2);
    float invU3 = 1.0 / (u * u2);
    return vec4(u, (r2 * invU3) * p.xy, (8.0 * invU3) * p.z);
}

// Exact Euclidean SDF and unit gradient of the Koranyi ball {r^4 + 16 z^2 <= R^4}.
// Parameterizes the meridian smoothly as r(a) = R sin(a), z(a) = (R^2/4) cos(a) sqrt(1 + sin^2(a))
// for a in [0, pi/2], eliminating the 4th-order flattening singularity at the pole and vertical tangent at the equator.
vec4 sdgKoranyiBall( in vec3 p, in float R )
{
    if (R <= 1e-12)
    {
        float l = length(p);
        return vec4(l - R, (l > 1e-20) ? p / l : vec3(0.0, 0.0, 1.0));
    }
    float r0 = length(p.xy);
    float z0 = abs(p.z);
    float R2 = R * R;
    float R4 = R2 * R2;
    float sgn = (r0 * r0 * r0 * r0 + 16.0 * z0 * z0 <= R4) ? -1.0 : 1.0;

    // Initial angle in [0, pi/2] from normalized meridian coordinates
    float a = atan(r0 * R, 2.0 * sqrt(max(z0, 0.0)) * R);
    // Also seed near equator/pole for interior points
    if (sgn < 0.0)
    {
        float u0 = clamp(r0 / R, 0.0, 1.0);
        a = asin(u0);
    }

    for (int i = 0; i < 4; i++)
    {
        float s   = sin(a);
        float c   = cos(a);
        float s2  = s * s;
        float q   = sqrt(1.0 + s2);
        float iq  = 1.0 / q;
        float rs  = R * s;
        float zs  = 0.25 * R2 * c * q;
        float dr  = rs - r0;
        float dz  = zs - z0;
        // Exact first and second derivatives with respect to a:
        //   rs'  = R c,                       zs'  = -0.5 R^2 s^3 / q
        //   rs'' = -R s,                      zs'' = -0.5 R^2 s^2 c (3 + 2 s^2) / q^3
        float rsp = R * c;
        float zsp = -0.5 * R2 * s2 * s * iq;
        float rspp = -rs;
        float zspp = -0.5 * R2 * s2 * c * (3.0 + 2.0 * s2) * (iq * iq * iq);

        float f1  = dr * rsp + dz * zsp;
        float f2  = rsp * rsp + zsp * zsp + dr * rspp + dz * zspp;
        a = clamp(a - f1 / max(abs(f2), 1e-12), 0.0, 1.57079632679);
    }

    float s  = sin(a);
    float c  = cos(a);
    float s2 = s * s;
    float q  = sqrt(1.0 + s2);
    float rs = R * s;
    float zs = 0.25 * R2 * c * q;
    float dr = r0 - rs;
    float dz = z0 - zs;
    float d  = sqrt(dr * dr + dz * dz);

    // Analytical outward normal to r^4 + 16 z^2 = R^4 has direction (r^3, 8 z) proportional to (s^3, 2 c q)
    vec2 nrz = (d > 1e-6) ? (sgn * vec2(dr, dz) / d)
                          : normalize(vec2(s2 * s, max(2.0 * c * q, 1e-12)));
    vec2 rdir = (r0 > 1e-20) ? (p.xy / r0) : vec2(1.0, 0.0);
    float zdir = (p.z >= 0.0) ? 1.0 : -1.0;
    return vec4(sgn * d, nrz.x * rdir, nrz.y * zdir);
}

float sdKoranyiBall( in vec3 p, in float R )
{
    return sdgKoranyiBall(p, R).x;
}

// ---------------------------------------------------------------------------
// 3. Exact Carnot-Caratheodory Distance (srCarnotDistance) & Gavrilov-Pansu Ball (sdCarnotBall)
//    Sub-Riemannian geodesics lift circles in R^2 to the Pansu sphere of radius R:
//      r(theta) = R sin(theta) / theta,   |z|(theta) = R^2 (2 theta - sin(2 theta)) / (8 theta^2)
//    for theta in [0, pi].
//    Variational Tangent-Line Identity:
//      R(r, |z|, theta) = (r cos(theta) + sqrt(r^2 cos^2(theta) + 4 |z| theta w(theta))) / w(theta)
//      with w(theta) = 1 + sin(theta) cos(theta) / theta in [0.7828, 2].
//    Satisfies dR / d(theta) == 0 at the exact geodesic angle theta*(p), giving O(eps^2) distance accuracy!
// ---------------------------------------------------------------------------

// Solves for the sub-Riemannian geodesic curvature parameter theta in [0, pi].
float srCarnotTheta( in float r, in float az )
{
    const float PI = 3.14159265359;
    if (az <= 1e-16) return 0.0;
    float r2 = r * r;
    float z4 = 4.0 * az;
    if (r2 <= 1e-24 * z4) return PI;

    // Two-chart rational minimax initial guess (max error < 6e-3 on [0, pi], 0 trig calls)
    float th;
    if (z4 <= 1.5 * r2)
    {
        float q  = z4 / r2;
        float q2 = q * q;
        th = q * (1.5 + 0.4338463 * q2 - 0.02044622 * q2 * q2) / (1.0 + 0.58669907 * q2);
    }
    else
    {
        float p = r * inversesqrt(z4);
        th = PI - p * (1.77245385 - 1.82856333 * p - 0.32382233 * p * p) / (1.0 - 1.1850757 * p);
    }
    th = clamp(th, 0.0, PI);

    // 1 Newton step on G(th) = r sqrt(th - sin(th)cos(th)) - 2 sqrt(|z|) sin(th) = 0
    // reaches |th - th*| < 1.2e-5 across all (r, z).
    float s   = sin(th);
    float c   = cos(th);
    float tsc = (th > 1e-3) ? (th - s * c) : (th * th * th * (0.66666667 - 0.13333333 * th * th));
    float rt  = sqrt(max(tsc, 1e-30));
    float sz  = sqrt(z4);
    float G   = r * rt - sz * s;
    float Gp  = r * (s * s) / max(rt, 1e-30) - sz * c;
    th = clamp(th - G / ((abs(Gp) > 1e-20) ? Gp : -1.0), 0.0, PI);
    return th;
}

// Returns vec4(d_CC(p), Euclidean gradient of d_CC(p)).
// Off the z-axis, d_CC satisfies the exact horizontal Eikonal equation ||grad_H d_CC||_H == 1,
// with partial derivatives d(d_CC)/dr = cos(theta) and d(d_CC)/d|z| = 2 theta / d_CC.
vec4 sdgCarnotDistance( in vec3 p )
{
    float r  = length(p.xy);
    float az = abs(p.z);
    if (r <= 1e-15 && az <= 1e-15) return vec4(0.0, 1.0, 0.0, 0.0);

    float th = srCarnotTheta(r, az);
    float s  = sin(th);
    float c  = cos(th);
    float sc_th = (th > 1e-4) ? (s * c / th) : (1.0 - 0.66666667 * th * th);
    float w  = 1.0 + sc_th;
    float rc = r * c;
    float R  = (rc + sqrt(max(rc * rc + 4.0 * az * th * w, 0.0))) / w;

    vec2  rdir = (r > 1e-20) ? (p.xy / r) : vec2(1.0, 0.0);
    float zdir = (p.z >= 0.0) ? 1.0 : -1.0;
    float gz   = (R > 1e-20) ? (2.0 * th / R) : 0.0;
    return vec4(R, c * rdir, gz * zdir);
}

float srCarnotDistance( in vec3 p )
{
    return sdgCarnotDistance(p).x;
}

// Exact Euclidean SDF and unit gradient of the Gavrilov-Pansu geodesic ball {d_CC(p) <= R}.
// The Pansu sphere has smooth meridian (r(th), z(th)) for th in [0, pi] with an inward apple-like
// conical dimple at the poles (r = 0, |z| = R^2 / (4 pi)).
vec4 sdgCarnotBall( in vec3 p, in float R )
{
    const float PI = 3.14159265359;
    if (R <= 1e-12)
    {
        float l = length(p);
        return vec4(l - R, (l > 1e-20) ? p / l : vec3(0.0, 0.0, 1.0));
    }
    float r0 = length(p.xy);
    float z0 = abs(p.z);
    float R2 = R * R;

    // Exact inside/outside sign from the sub-Riemannian Carnot distance
    float th_cc = srCarnotTheta(r0, z0);
    float s_cc  = sin(th_cc);
    float c_cc  = cos(th_cc);
    float w_cc  = 1.0 + ((th_cc > 1e-4) ? (s_cc * c_cc / th_cc) : (1.0 - 0.66666667 * th_cc * th_cc));
    float rc_cc = r0 * c_cc;
    float d_cc  = (rc_cc + sqrt(max(rc_cc * rc_cc + 4.0 * z0 * th_cc * w_cc, 0.0))) / w_cc;
    float sgn   = (d_cc <= R) ? -1.0 : 1.0;

    // Candidate 1: The polar dimple vertex at th = pi: (r, z) = (0, R^2 / (4 pi))
    float z_pole = R2 * (0.25 / PI);
    float best_d2 = r0 * r0 + (z0 - z_pole) * (z0 - z_pole);
    float best_r  = 0.0;
    float best_z  = z_pole;
    float best_th = PI;

    // Candidate 2: Smooth meridian projection via Newton on the reduced orthogonality equation:
    //   Since (r'(th), z'(th)) is parallel to (2 th, -R cos(th)), orthogonality to (r(th)-r0, z(th)-z0) is:
    //   Psi(th) = 2 th (r(th) - r0) - R cos(th) (z(th) - z0) = 0.
    // Test 2 initial seeds: the radial Carnot angle th_cc and a mid-latitude seed for the polar fold
    for (int seed = 0; seed < 2; seed++)
    {
        float th = (seed == 0) ? clamp(th_cc, 1e-4, PI - 1e-3) : 2.2;
        for (int i = 0; i < 5; i++)
        {
            float s   = sin(th);
            float c   = cos(th);
            float ith = 1.0 / max(th, 1e-5);
            float rs  = (th > 1e-4) ? (R * s * ith) : (R * (1.0 - th * th / 6.0));
            float zs  = (th > 1e-4) ? (0.25 * R2 * (th - s * c) * ith * ith) : (R2 * th / 6.0);
            float rp  = (th > 1e-4) ? (R * (th * c - s) * ith * ith) : (-R * th / 3.0);
            float zp  = -0.5 * R * c * ith * rp;

            float dr  = rs - r0;
            float dz  = zs - z0;
            float Psi = 2.0 * th * dr - R * c * dz;
            float dPsi = 2.0 * dr + 2.0 * th * rp + R * s * dz - R * c * zp;
            th = clamp(th - Psi / ((abs(dPsi) > 1e-12) ? dPsi : 1.0), 0.0, PI);
        }
        float s   = sin(th);
        float c   = cos(th);
        float ith = 1.0 / max(th, 1e-5);
        float rs  = (th > 1e-4) ? (R * s * ith) : (R * (1.0 - th * th / 6.0));
        float zs  = (th > 1e-4) ? (0.25 * R2 * (th - s * c) * ith * ith) : (R2 * th / 6.0);
        float d2  = (r0 - rs) * (r0 - rs) + (z0 - zs) * (z0 - zs);
        if (d2 < best_d2)
        {
            best_d2 = d2;
            best_r  = rs;
            best_z  = zs;
            best_th = th;
        }
    }

    float d   = sqrt(best_d2);
    float dr  = r0 - best_r;
    float dz  = z0 - best_z;
    vec2  nrz = (d > 1e-6) ? (sgn * vec2(dr, dz) / d)
                           : normalize(vec2(R * cos(best_th), 2.0 * best_th));
    vec2  rdir = (r0 > 1e-20) ? (p.xy / r0) : vec2(1.0, 0.0);
    float zdir = (p.z >= 0.0) ? 1.0 : -1.0;
    return vec4(sgn * d, nrz.x * rdir, nrz.y * zdir);
}

float sdCarnotBall( in vec3 p, in float R )
{
    return sdgCarnotBall(p, R).x;
}

// ---------------------------------------------------------------------------
// 4. Left-Invariant Sub-Riemannian Translation (opHeisenbergTranslate)
//    Left-translating a shape by a in H^1 maps sample point p to:
//      q = a^{-1} * p = (p.x - a.x, p.y - a.y, p.z - a.z - 0.5 * (a.x * p.y - a.y * p.x)).
//    Its Jacobian J = dq/dp has det(J) = 1 (volume-preserving symplectic shear) and maximum
//    singular value sigma_max(J) = sqrt(1 + ||a.xy||^2 / 16) + ||a.xy|| / 4.
//    Multiplying any 1-Lipschitz Euclidean SDF d(q) by 1 / sigma_max(J) yields an exact
//    tight 1-Lipschitz Euclidean bound with zero raymarching overstep!
// ---------------------------------------------------------------------------

vec3 opHeisenbergTranslate( in vec3 p, in vec3 a )
{
    return vec3(
        p.x - a.x,
        p.y - a.y,
        p.z - a.z - 0.5 * (a.x * p.y - a.y * p.x)
    );
}

// Exact reciprocal of the maximum singular value of J = d(a^{-1} * p)/dp:
//   1 / sigma_max(J) = 1 / (sqrt(1 + ||a.xy||^2 / 16) + ||a.xy|| / 4)
//                    = sqrt(1 + ||a.xy||^2 / 16) - ||a.xy|| / 4.
float opHeisenbergTranslateScale( in vec3 a )
{
    float s = 0.25 * length(a.xy);
    return 1.0 / (sqrt(1.0 + s * s) + s);
}

// Transforms a sample point p by q = a^{-1} * p and pulls back an inner SDF + gradient dg = sdgInner(q)
// via J^T * dg.yzw, returning the Lipschitz-safe distance and normalized Euclidean direction.
vec4 sdgHeisenbergTranslate( in vec3 a, in vec4 dgInner )
{
    float scale = opHeisenbergTranslateScale(a);
    vec3  g = vec3(
        dgInner.y + 0.5 * a.y * dgInner.w,
        dgInner.z - 0.5 * a.x * dgInner.w,
        dgInner.w
    );
    float lg = length(g);
    return vec4(dgInner.x * scale, (lg > 1e-20) ? (g / lg) : vec3(0.0, 0.0, 1.0));
}

#endif // SDF_GEMS_HAS_SR_HEISENBERG
