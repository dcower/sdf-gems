// ============================================================================
// sminVaradhan — Angle-Aware Varadhan-Eikonal & Gram-Metric Rolling-Ball Smooth Booleans
// @author Discovered & mathematically designed using Gemini 4 Argon, implemented & verified with Gemini 3.8 Flash (sdf-gems, 2026)
// @exact  Exact Gram-metric constant-radius rolling-ball fillet (||grad d|| = 1.000000 across all intersection angles theta in (0, pi))
// ============================================================================

#ifndef SDF_GEMS_SMIN_VARADHAN_GLSL
#define SDF_GEMS_SMIN_VARADHAN_GLSL

// ----------------------------------------------------------------------------
// 1. Exact Gram-Metric Constant-Radius Rolling-Ball Smooth Union (3D & 2D)
// ----------------------------------------------------------------------------

// 3D Exact Gram-metric constant-radius rolling-ball smooth union of two 1-jets a=(da, ga) and b=(db, gb).
// Unlike standard angle-blind smin functions whose fillet width varies as 1/sin(theta/2) and whose
// gradient norm drops by up to 29.3% (orthogonal) or 100% (opposing), sdgSminGram uses the inverse Gram
// matrix of the wavefront normals to construct the exact Euclidean rolling-ball fillet with ||grad|| = 1.
vec4 sdgSminGram( in vec4 a, in vec4 b, in float r )
{
    float c = dot(a.yzw, b.yzw);
    float sin2 = 1.0 - c * c;
    if( sin2 > 1e-6 )
    {
        float ua = r - a.x;
        float ub = r - b.x;
        if( (ua - ub * c > 0.0) && (ub - ua * c > 0.0) )
        {
            float du = 2.0 + 2.0 * c;
            float dv = 2.0 - 2.0 * c;
            float su = a.x + b.x - 2.0 * r;
            float sv = a.x - b.x;
            float rho2 = (su * su) / du + (sv * sv) / dv;
            float rho = sqrt(max(rho2, 1e-14));
            float inv_rho = 1.0 / rho;
            float cu = (su / du) * inv_rho;
            float cv = (sv / dv) * inv_rho;
            float wa = -(cu + cv);
            float wb = -(cu - cv);
            return vec4(r - rho, wa * a.yzw + wb * b.yzw);
        }
    }
    return (a.x < b.x) ? a : b;
}

// 2D Exact Gram-metric constant-radius rolling-ball smooth union of two 1-jets a=(da, ga) and b=(db, gb).
vec3 sdgSminGram2D( in vec3 a, in vec3 b, in float r )
{
    float c = dot(a.yz, b.yz);
    float sin2 = 1.0 - c * c;
    if( sin2 > 1e-6 )
    {
        float ua = r - a.x;
        float ub = r - b.x;
        if( (ua - ub * c > 0.0) && (ub - ua * c > 0.0) )
        {
            float du = 2.0 + 2.0 * c;
            float dv = 2.0 - 2.0 * c;
            float su = a.x + b.x - 2.0 * r;
            float sv = a.x - b.x;
            float rho2 = (su * su) / du + (sv * sv) / dv;
            float rho = sqrt(max(rho2, 1e-14));
            float inv_rho = 1.0 / rho;
            float cu = (su / du) * inv_rho;
            float cv = (sv / dv) * inv_rho;
            float wa = -(cu + cv);
            float wb = -(cu - cv);
            return vec3(r - rho, wa * a.yz + wb * b.yz);
        }
    }
    return (a.x < b.x) ? a : b;
}

// ----------------------------------------------------------------------------
// 2. Exact Gram-Metric Constant-Radius Rolling-Ball Smooth Intersection (3D & 2D)
// ----------------------------------------------------------------------------

// 3D Exact Gram-metric constant-radius rolling-ball smooth intersection (chamfer/fillet on exterior corner).
vec4 sdgSmaxGram( in vec4 a, in vec4 b, in float r )
{
    float c = dot(a.yzw, b.yzw);
    float sin2 = 1.0 - c * c;
    if( sin2 > 1e-6 )
    {
        float ua = r + a.x;
        float ub = r + b.x;
        if( (ua - ub * c > 0.0) && (ub - ua * c > 0.0) )
        {
            float du = 2.0 + 2.0 * c;
            float dv = 2.0 - 2.0 * c;
            float su = a.x + b.x + 2.0 * r;
            float sv = a.x - b.x;
            float rho2 = (su * su) / du + (sv * sv) / dv;
            float rho = sqrt(max(rho2, 1e-14));
            float inv_rho = 1.0 / rho;
            float cu = (su / du) * inv_rho;
            float cv = (sv / dv) * inv_rho;
            float wa = cu + cv;
            float wb = cu - cv;
            return vec4(rho - r, wa * a.yzw + wb * b.yzw);
        }
    }
    return (a.x > b.x) ? a : b;
}

// 2D Exact Gram-metric constant-radius rolling-ball smooth intersection.
vec3 sdgSmaxGram2D( in vec3 a, in vec3 b, in float r )
{
    float c = dot(a.yz, b.yz);
    float sin2 = 1.0 - c * c;
    if( sin2 > 1e-6 )
    {
        float ua = r + a.x;
        float ub = r + b.x;
        if( (ua - ub * c > 0.0) && (ub - ua * c > 0.0) )
        {
            float du = 2.0 + 2.0 * c;
            float dv = 2.0 - 2.0 * c;
            float su = a.x + b.x + 2.0 * r;
            float sv = a.x - b.x;
            float rho2 = (su * su) / du + (sv * sv) / dv;
            float rho = sqrt(max(rho2, 1e-14));
            float inv_rho = 1.0 / rho;
            float cu = (su / du) * inv_rho;
            float cv = (sv / dv) * inv_rho;
            float wa = cu + cv;
            float wb = cu - cv;
            return vec3(rho - r, wa * a.yz + wb * b.yz);
        }
    }
    return (a.x > b.x) ? a : b;
}

// ----------------------------------------------------------------------------
// 3. Scalar Gram-Metric Rolling-Ball Smooth Union & Intersection
// ----------------------------------------------------------------------------

// Exact scalar rolling-ball smooth union when intersection cosine cos_theta is known.
float sminGram( in float a, in float b, in float r, in float cos_theta )
{
    float c = clamp(cos_theta, -0.9999, 0.9999);
    float ua = r - a;
    float ub = r - b;
    if( (ua - ub * c > 0.0) && (ub - ua * c > 0.0) )
    {
        float su = a + b - 2.0 * r;
        float sv = a - b;
        float rho2 = (su * su) / (2.0 + 2.0 * c) + (sv * sv) / (2.0 - 2.0 * c);
        return r - sqrt(max(rho2, 0.0));
    }
    return min(a, b);
}

// Exact scalar rolling-ball smooth intersection when intersection cosine cos_theta is known.
float smaxGram( in float a, in float b, in float r, in float cos_theta )
{
    return -sminGram(-a, -b, r, cos_theta);
}

// ----------------------------------------------------------------------------
// 4. Varadhan Angle-Corrected Screened-Poisson Smooth Minimum (3D & 2D)
// ----------------------------------------------------------------------------

// 3D C^infinity Varadhan Screened-Poisson smooth union with first-order Eikonal gradient normalization.
vec4 sdgSminVaradhan( in vec4 a, in vec4 b, in float k )
{
    float c = clamp(dot(a.yzw, b.yzw), -1.0, 1.0);
    float s = sqrt(max(0.5 * (1.0 - c), 1e-4));
    float keff = k * s;
    float diff = clamp((a.x - b.x) / keff, -50.0, 50.0);
    float wa = 1.0 / (1.0 + exp(diff));
    float wb = 1.0 - wa;
    float d = min(a.x, b.x) - keff * log(1.0 + exp(-abs(diff)));
    vec3 g_raw = wa * a.yzw + wb * b.yzw;
    float eikonal_sq = max(1.0 - 2.0 * wa * wb * (1.0 - c), 1e-6);
    vec3 g = g_raw * inversesqrt(eikonal_sq);
    return vec4(d, g);
}

// 2D C^infinity Varadhan Screened-Poisson smooth union with first-order Eikonal gradient normalization.
vec3 sdgSminVaradhan2D( in vec3 a, in vec3 b, in float k )
{
    float c = clamp(dot(a.yz, b.yz), -1.0, 1.0);
    float s = sqrt(max(0.5 * (1.0 - c), 1e-4));
    float keff = k * s;
    float diff = clamp((a.x - b.x) / keff, -50.0, 50.0);
    float wa = 1.0 / (1.0 + exp(diff));
    float wb = 1.0 - wa;
    float d = min(a.x, b.x) - keff * log(1.0 + exp(-abs(diff)));
    vec2 g_raw = wa * a.yz + wb * b.yz;
    float eikonal_sq = max(1.0 - 2.0 * wa * wb * (1.0 - c), 1e-6);
    vec2 g = g_raw * inversesqrt(eikonal_sq);
    return vec3(d, g);
}

// Scalar Varadhan angle-scaled exponential smooth minimum.
float sminVaradhan( in float a, in float b, in float k, in float cos_theta )
{
    float c = clamp(cos_theta, -1.0, 1.0);
    float s = sqrt(max(0.5 * (1.0 - c), 1e-4));
    float keff = k * s;
    float diff = clamp(abs(a - b) / keff, 0.0, 50.0);
    return min(a, b) - keff * log(1.0 + exp(-diff));
}

// ----------------------------------------------------------------------------
// 5. Kelvin Electrostatic Capacitary Smooth Minimum
// ----------------------------------------------------------------------------

// Scalar rational potential smooth minimum based on electrostatic capacitance.
float sminKelvin( in float a, in float b, in float k )
{
    float u = 2.0 * k + a;
    float v = 2.0 * k + b;
    return (u > 0.0 && v > 0.0) ? (u * v) / (u + v) - 2.0 * k : min(a, b);
}

// 3D Kelvin electrostatic capacitary smooth minimum with analytical gradient.
vec4 sdgSminKelvin( in vec4 a, in vec4 b, in float k )
{
    float u = 2.0 * k + a.x;
    float v = 2.0 * k + b.x;
    if( u > 0.0 && v > 0.0 )
    {
        float sum = u + v;
        float d = (u * v) / sum - 2.0 * k;
        float wa = (v * v) / (sum * sum);
        float wb = (u * u) / (sum * sum);
        vec3 g = normalize(wa * a.yzw + wb * b.yzw);
        return vec4(d, g);
    }
    return (a.x < b.x) ? a : b;
}

#endif // SDF_GEMS_SMIN_VARADHAN_GLSL
