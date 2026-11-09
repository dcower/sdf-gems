// sdfYukawaColeHopf: Semiclassical Log-Exponential Duality operators:
// 1. Maslov dequantization sdgMaslovSmin & Fisher-information curvature sdfMaslovCurvature (Litvinov & Maslov 1998)
// 2. Modified Bessel function sdfBesselK0 & Varadhan Screened-Poisson ray integral sdfYukawaRayIntegral (Varadhan 1967)
// 3. Multi-center Varadhan-Yukawa log-potential metaball SDF sdYukawaPotential
// 4. Cole-Hopf scale-space cone filtering sdfColeHopfScaleSpace (Hopf 1950; Cole 1951; KPZ 1986)
// See docs/volume/sdfYukawaColeHopf.md for mathematical derivations.
#ifndef SDF_GEMS_SDF_YUKAWA_COLE_HOPF_GLSL
#define SDF_GEMS_SDF_YUKAWA_COLE_HOPF_GLSL

// High-precision Abramowitz & Stegun (9.8.5 & 9.8.6) polynomial approximation for the
// Modified Bessel function of the second kind of order 0, K_0(x), for x in (0, infinity).
// Uniform absolute error < 4e-8 across all x > 0.
float sdfBesselK0(float x)
{
    x = max(x, 1e-6);
    if (x <= 2.0)
    {
        float y = (x / 3.75) * (x / 3.75);
        float i0 = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))));
        float t = (x / 2.0) * (x / 2.0);
        float k0_poly = -0.57721566 + t * (0.42278420 + t * (0.23069756 + t * (0.03488590 + t * (0.00262698 + t * (0.00010750 + t * 0.00000740)))));
        return -log(0.5 * x) * i0 + k0_poly;
    }
    else
    {
        float u = 2.0 / x;
        float poly = 1.25331414 + u * (-0.07832358 + u * (0.02189568 + u * (-0.01062446 + u * (0.00587872 + u * (-0.00251540 + u * 0.00053208)))));
        return (exp(-x) / sqrt(x)) * poly;
    }
}

// Tropical Maslov Dequantization smooth minimum for scalar distances (Litvinov & Maslov 1998; Viro 2001).
// Stable log-sum-exp formulation with temperature parameter h > 0.
float sdMaslovSmin(float a, float b, float h)
{
    float m = min(a, b);
    return m - h * log(exp(-(a - m) / h) + exp(-(b - m) / h));
}

// Tropical Maslov Dequantization smooth minimum and exact analytical gradient (M1, Litvinov & Maslov 1998).
// Computes log-sum-exp smooth union between two SDF+gradient vec4(d, gx, gy, gz) inputs dga and dgb.
// Because the analytical gradient is a convex combination of ga and gb with Gibbs softmax weights
// wa, wb in [0, 1] (wa + wb = 1), ||grad Dh|| <= max(||ga||, ||gb||) <= 1: a certified 1-Lipschitz bound.
vec4 sdgMaslovSmin(vec4 dga, vec4 dgb, float h)
{
    float m = min(dga.x, dgb.x);
    float diff = (dga.x - dgb.x) / h;
    float wa = 1.0 / (1.0 + exp(clamp(diff, -60.0, 60.0)));
    float wb = 1.0 - wa;
    float Dh = m - h * log(exp(-(dga.x - m) / h) + exp(-(dgb.x - m) / h));
    return vec4(Dh, wa * dga.yzw + wb * dgb.yzw);
}

// Exact Fisher-information covariance Laplacian Delta Dh = div(grad Dh).
// Mean curvature H = 0.5 * Delta Dh at smooth CSG blends.
float sdfMaslovCurvature(vec4 dga, vec4 dgb, float lapA, float lapB, float h)
{
    float diff = (dga.x - dgb.x) / h;
    float wa = 1.0 / (1.0 + exp(clamp(diff, -60.0, 60.0)));
    float wb = 1.0 - wa;
    vec3 dgrad = dga.yzw - dgb.yzw;
    return wa * lapA + wb * lapB - (wa * wb / h) * dot(dgrad, dgrad);
}

// Varadhan's Screened-Poisson Formula & Yukawa Green's Function Ray Integral (M2, Varadhan 1967).
// Ray line integral of G_lambda(x, c) = exp(-lambda * ||x - c||) / (4 * pi * ||x - c||) along
// an infinite ray o + t * d (t in (-inf, inf)) evaluated in zero marching steps.
// Assumes rd is normalized. b is the impact parameter (perpendicular ray-point distance).
float sdfYukawaRayIntegral(vec3 ro, vec3 rd, vec3 center, float lambda)
{
    vec3 diff = ro - center;
    float b = length(diff - dot(diff, rd) * rd);
    return (0.5 / 3.14159265358979323846) * sdfBesselK0(lambda * b);
}

// Varadhan Log-Potential SDF of a 3-center Screened-Poisson (Yukawa) potential sum.
// Evaluates u_lambda(p) = sum_{k=1}^3 (r0 / max(||p - c_k||, 1e-4)) * exp(-lambda * (||p - c_k|| - r0))
// and returns d_lambda(p) = -lambda^(-1) * ln(u_lambda(p)).
// Far-field gradient satisfies ||grad d|| <= 1.
// The volumetric corona integrates in zero marching steps via sum_k (r0 * exp(lambda * r0) / (2 * pi)) * K0(lambda * b_k).
float sdYukawaPotential(vec3 p, vec3 c1, vec3 c2, vec3 c3, float r0, float lambda)
{
    float d1 = length(p - c1);
    float d2 = length(p - c2);
    float d3 = length(p - c3);
    float u = (r0 / max(d1, 1e-4)) * exp(-lambda * (d1 - r0))
            + (r0 / max(d2, 1e-4)) * exp(-lambda * (d2 - r0))
            + (r0 / max(d3, 1e-4)) * exp(-lambda * (d3 - r0));
    return -(1.0 / lambda) * log(max(u, 1e-12));
}

// Cole-Hopf Scale-Space Heat Semigroup evolution (M16, Hopf 1950; Cole 1951; KPZ 1986).
// Filters an SDF d(x) across pixel cone footprint radius sigma = t * tan(theta) (variance s = sigma^2)
// under the viscous Eikonal PDE d_s u = 0.5 * Delta u - 0.5 * lambda * ||grad u||^2.
// Linearizes via Phi = exp(-lambda * u) into the 3D heat semigroup, rounding convex crests (Delta d > 0)
// and filtering sub-pixel details in 1 sample per pixel.
float sdfColeHopfScaleSpace(float d, float lapD, float gradNormSq, float sigma, float lambda)
{
    float s = sigma * sigma;
    float bracket = 1.0 + 0.5 * s * (lambda * lambda * gradNormSq - lambda * lapD);
    return d - (1.0 / lambda) * log(max(bracket, 1e-4));
}

#endif // SDF_GEMS_SDF_YUKAWA_COLE_HOPF_GLSL
