// ============================================================================
// rmEigenRepulsion — 3x3 Hermitian Eigenvalue Repulsion, Secular Smooth Booleans,
//                    Hellmann-Feynman Gradients & Weyl Diabolical Cones
//   rmEigenStates3, rmEigenWeights3, opEigenUnion3, opEigenIntersection3,
//   sdgEigenUnion3, sdgEigenTripod3D, sdEigenTripod3D,
//   sdEigenExcitedShell3D, sdVandermondeTube3D
// @author Gemini 4 Argon (sdf-gems, 2026)
// @exact  Every eigenvalue sheet lambda_m(p) (m = 1, 2, 3) and the half-gap
//         0.5 * (lambda_2(p) - lambda_1(p)) are globally 1-Lipschitz conservative
//         SDF bounds across all of R^3 for ANY 1-Lipschitz inputs d_1, d_2, d_3
//         and ANY constant couplings k = (k_12, k_13, k_23) (by Weyl's inequality
//         and the Hellmann-Feynman convex-combination theorem).
// ============================================================================

#ifndef SDF_GEMS_HAS_RM_EIGEN_REPULSION
#define SDF_GEMS_HAS_RM_EIGEN_REPULSION

// ----------------------------------------------------------------------------
// 1. Branchless 3x3 Real Symmetric Eigenvalue Solver (5D Grassmannian 2-Blade)
//    Computes all three ordered eigenvalues (lambda_1 <= lambda_2 <= lambda_3) of
//        [ d.x   k.x   k.y ]
//    H = [ k.x   d.y   k.z ]
//        [ k.y   k.z   d.z ]
//    By shifting B = H - mu * I (mu = tr(H)/3) and mapping B and Y = B^2 - (tr(B^2)/3)*I
//    into the 5D Euclidean space of traceless symmetric matrices (Cartan R^2 + root R^3),
//    Viete's angle 3*theta is the exact Euclidean angle between B and Y:
//        3*theta = atan( ||b ^ y||_2, <b, y> ).
//    Evaluating ||b ^ y||_2 as a sum of squares of 2x2 minors eliminates the O(sqrt(eps))
//    catastrophic cancellation of textbook 4*p^6 - det(B)^2 solvers near Weyl degeneracies
//    and in the sharp k -> 0 boolean limit.
// ----------------------------------------------------------------------------
vec3 rmEigenStates3( in vec3 d, in vec3 k )
{
    const float s3 = 1.7320508075688772;
    float mu  = (d.x + d.y + d.z) * (1.0 / 3.0);
    float d12 = d.x - d.y;
    float d13 = d.x - d.z;
    float d23 = d.y - d.z;

    float a1 = (d12 + d13) * (1.0 / 3.0);
    float a2 = (d23 - d12) * (1.0 / 3.0);
    float a3 = -(d13 + d23) * (1.0 / 3.0);

    // 5D traceless coordinates of B = (b01 in Cartan R^2, b234 in root R^3), ||b|| = 3*p
    vec2 b01  = vec2(0.5 * s3 * d12, 0.5 * (d13 + d23));
    vec3 b234 = s3 * k;

    // 5D traceless coordinates of Y = B^2 - (tr(B^2)/3)*I
    float kd0 = 0.5 * s3 * (k.y - k.z) * (k.y + k.z);
    float kd1 = 0.5 * (2.0 * k.x * k.x - k.y * k.y - k.z * k.z);
    vec2  y01 = vec2(-0.5 * s3 * a3 * d12 + kd0,
                      0.5 * ((a1 - a3) * (a1 + a3) + (a2 - a3) * (a2 + a3)) + kd1);
    vec3  y234 = s3 * vec3(-a3 * k.x + k.y * k.z,
                           -a2 * k.y + k.x * k.z,
                           -a1 * k.z + k.x * k.y);

    // Factored Cartan-root 2-blade norm ||b ^ y||_2 = (3*sqrt(3)/4) * VandermondeDiscriminant
    float w01   = -0.5 * s3 * d12 * d13 * d23 + b01.x * kd1 - b01.y * kd0;
    vec3  w0    = b01.x * y234 - y01.x * b234;
    vec3  w1    = b01.y * y234 - y01.y * b234;
    vec3  w234  = cross(b234, y234);
    float wedge = sqrt(w01 * w01 + dot(w0, w0) + dot(w1, w1) + dot(w234, w234));

    float dot_by = dot(b01, y01) + dot(b234, y234);
    float p      = sqrt(dot(b01, b01) + dot(b234, b234)) * (1.0 / 3.0);
    float theta  = atan(wedge, dot_by + 1e-38) * (1.0 / 3.0);
    float c      = cos(theta);
    float s      = sin(theta);

    return vec3(mu - p * (c + s3 * s),
                mu - p * (c - s3 * s),
                mu + 2.0 * p * c);
}

vec3 rmEigenStates3( in vec3 d, in float k )
{
    return rmEigenStates3(d, vec3(k));
}

// ----------------------------------------------------------------------------
// 2. Eigenvalue-Eigenvector Identity (Adjugate Principal Minors, 0 sqrt)
//    Returns the Hellmann-Feynman probability weights w_i = v_{m,i}^2 >= 0
//    (with w.x + w.y + w.z = 1) for eigenvalue lambda_m, using the identity
//        adj(H - lambda_m I)_{ii} = (lambda_m - lambda_a)(lambda_m - lambda_b) v_{m,i}^2.
//    When lambda_m approaches a double eigenvalue (all principal 2x2 minors -> 0),
//    smoothly falls back to 0.5 * (1 - w_other) via the trace identity P_1 + P_2 + P_3 = I.
// ----------------------------------------------------------------------------
vec3 rmEigenWeights3( in vec3 d, in vec3 k, in float lambda_m, in float lambda_other )
{
    vec3  dl = d - vec3(lambda_m);
    vec3  k2 = k * k;
    vec3  m  = abs(vec3(dl.y * dl.z - k2.z,
                        dl.x * dl.z - k2.y,
                        dl.x * dl.y - k2.x));
    float sm = m.x + m.y + m.z;

    vec3  dlo = d - vec3(lambda_other);
    vec3  mo  = abs(vec3(dlo.y * dlo.z - k2.z,
                         dlo.x * dlo.z - k2.y,
                         dlo.x * dlo.y - k2.x));
    float so  = mo.x + mo.y + mo.z;

    vec3 w_other = (so > 1e-20) ? (mo / so) : vec3(1.0 / 3.0);
    vec3 w_deg   = max(0.5 * (vec3(1.0) - w_other), vec3(0.0));
    return (sm > 1e-12 * (so + 1e-20)) ? (m / sm) : w_deg;
}

// ----------------------------------------------------------------------------
// 3. S_3-Symmetric Secular Smooth Booleans & Hellmann-Feynman Gradients
//    - Unfrustrated Perron-Frobenius coupling (k_ij = -|k| < 0, k_12*k_13*k_23 < 0):
//      the ground state lambda_1 is strictly simple everywhere in R^3 (Delta_12 >= 2|k|),
//      yielding a C^infinity S_3-symmetric smooth union with zero triple-junction crease.
//    - Top state lambda_3 with frustrated coupling (+|k|) equals -lambda_1(-d, -|k|),
//      yielding the dual C^infinity S_3-symmetric smooth intersection.
// ----------------------------------------------------------------------------
float opEigenUnion3( in vec3 d, in vec3 k )
{
    return rmEigenStates3(d, k).x;
}

float opEigenUnion3( in vec3 d, in float k )
{
    return rmEigenStates3(d, vec3(-abs(k))).x;
}

float opEigenIntersection3( in vec3 d, in vec3 k )
{
    return rmEigenStates3(d, k).z;
}

float opEigenIntersection3( in vec3 d, in float k )
{
    return rmEigenStates3(d, vec3(abs(k))).z;
}

// Analytical value + unit normal via the Hellmann-Feynman theorem:
//     grad(lambda_1) = w.x * grad(d_1) + w.y * grad(d_2) + w.z * grad(d_3).
vec4 sdgEigenUnion3( in vec4 dg1, in vec4 dg2, in vec4 dg3, in vec3 k )
{
    vec3  d   = vec3(dg1.x, dg2.x, dg3.x);
    vec3  eig = rmEigenStates3(d, k);
    vec3  w   = rmEigenWeights3(d, k, eig.x, eig.z);
    vec3  g   = w.x * dg1.yzw + w.y * dg2.yzw + w.z * dg3.yzw;
    float l2  = dot(g, g);
    vec3  gn  = (l2 > 1e-24) ? (g * inversesqrt(l2)) : dg1.yzw;
    return vec4(eig.x, gn);
}

vec4 sdgEigenUnion3( in vec4 dg1, in vec4 dg2, in vec4 dg3, in float k )
{
    return sdgEigenUnion3(dg1, dg2, dg3, vec3(-abs(k)));
}

// ----------------------------------------------------------------------------
// 4. Concrete 3D Primitives: S_3 Tripod Cluster, Excited-State Shell, and
//    Vandermonde / Weyl Diabolical Cone Tube
// ----------------------------------------------------------------------------
vec4 sdgEigenSphereHelper( in vec3 p, in vec3 c, in float r )
{
    vec3  q  = p - c;
    float l2 = dot(q, q);
    if (l2 > 1e-24)
    {
        float invL = inversesqrt(l2);
        return vec4(l2 * invL - r, q * invL);
    }
    return vec4(-r, vec3(1.0, 0.0, 0.0));
}

// Three spheres of radius r at the vertices of an equilateral triangle of
// circumradius R in the xz-plane, blended by Perron-Frobenius eigenvalue repulsion.
vec4 sdgEigenTripod3D( in vec3 p, in float R, in float r, in float k )
{
    const float s3_2 = 0.8660254037844386;
    vec4 dg1 = sdgEigenSphereHelper(p, vec3( R,        0.0,  0.0     ), r);
    vec4 dg2 = sdgEigenSphereHelper(p, vec3(-0.5 * R,  0.0,  s3_2 * R), r);
    vec4 dg3 = sdgEigenSphereHelper(p, vec3(-0.5 * R,  0.0, -s3_2 * R), r);
    return sdgEigenUnion3(dg1, dg2, dg3, vec3(-abs(k)));
}

float sdEigenTripod3D( in vec3 p, in float R, in float r, in float k )
{
    return sdgEigenTripod3D(p, R, r, k).x;
}

// Excited-state saddle sheet lambda_2(p) (or thin shell |lambda_2(p)| - th when th > 0):
// globally 1-Lipschitz internal avoided-crossing chamber surface.
float sdEigenExcitedShell3D( in vec3 p, in float R, in float r, in float k, in float th )
{
    const float s3_2 = 0.8660254037844386;
    vec3 d = vec3(
        length(p - vec3( R,        0.0,  0.0     )) - r,
        length(p - vec3(-0.5 * R,  0.0,  s3_2 * R)) - r,
        length(p - vec3(-0.5 * R,  0.0, -s3_2 * R)) - r
    );
    float l2 = rmEigenStates3(d, vec3(-abs(k))).y;
    return (th > 0.0) ? (abs(l2) - th) : l2;
}

// Vandermonde / Weyl Diabolical Cone Tube:
// For frustrated coupling (+|k|), eigenvalue repulsion gaps out all three 2D pairwise
// bisector sheets (Delta_12 -> 2|k|) while topologically protecting the 1D conical
// degeneracy Delta_12 = 0 along the triple-equidistance Weyl locus d_1 = d_2 = d_3.
// Since lambda_1 and lambda_2 are each 1-Lipschitz, 0.5 * (lambda_2 - lambda_1) is
// globally 1-Lipschitz across R^3; clipping with lambda_1 preserves 1-Lipschitzness.
float sdVandermondeTube3D( in vec3 p, in float R, in float r, in float k, in float r_tube )
{
    const float s3_2 = 0.8660254037844386;
    vec3 d = vec3(
        length(p - vec3( R,        0.0,  0.0     )) - r,
        length(p - vec3(-0.5 * R,  0.0,  s3_2 * R)) - r,
        length(p - vec3(-0.5 * R,  0.0, -s3_2 * R)) - r
    );
    vec3  eig    = rmEigenStates3(d, vec3(abs(k)));
    float d_gap  = 0.5 * (eig.y - eig.x) - r_tube;
    return max(d_gap, eig.x);
}

#endif // SDF_GEMS_HAS_RM_EIGEN_REPULSION
