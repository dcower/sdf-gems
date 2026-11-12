#ifndef SDF_GEMS_HAS_BORROMEAN_RINGS
#define SDF_GEMS_HAS_BORROMEAN_RINGS

#include "../common/common.glsl"
#include "../primitives2d/sdEllipse2D.glsl"

// Forward declarations for single-translation-unit compilation in arbitrary order
float sdEllipse2D( in vec2 p, in vec2 ab );
vec3 sdgEllipse2D( in vec2 p, in vec2 ab );

// ============================================================================
// Exact 3D Elliptic Torus, Borromean Rings & 0-Loop Torus Knot
// @author Gemini 4 Argon (research & mathematical design) and Gemini 3.1 Pro (experiments & implementation)
// @exact  Yes (Elliptic Torus and Borromean Rings are exactly Euclidean)
// ============================================================================

// 1. Exact Euclidean 3D Elliptic Torus SDF
float sdEllipticTorus(vec3 p, vec2 ab, float r) {
    float d2 = sdEllipse2D(p.xy, ab);
    return sqrt(max(d2 * d2 + p.z * p.z, 0.0)) - r;
}

// 2. Exact Euclidean 3D Elliptic Torus SDF + Analytical Unit Gradient
vec4 sdgEllipticTorus(vec3 p, vec2 ab, float r) {
    vec3 d2 = sdgEllipse2D(p.xy, ab);
    float d = d2.x;
    float dist = sqrt(max(d * d + p.z * p.z, 0.0));
    vec3 grad = (dist > 1e-20) ? vec3(d2.yz * (d / dist), p.z / dist) : vec3(0.0, 0.0, 1.0);
    return vec4(dist - r, grad);
}

// 3. Exact Euclidean 3D Borromean Rings SDF
float sdBorromeanRings(vec3 p, vec2 ab, float r) {
    float d1 = sdEllipticTorus(p.xyz, ab, r);
    float d2 = sdEllipticTorus(p.yzx, ab, r);
    float d3 = sdEllipticTorus(p.zxy, ab, r);
    return min(min(d1, d2), d3);
}

// 4. Exact Euclidean 3D Borromean Rings SDF + Analytical Unit Gradient
vec4 sdgBorromeanRings(vec3 p, vec2 ab, float r) {
    vec4 d1 = sdgEllipticTorus(p.xyz, ab, r);
    vec4 d2 = sdgEllipticTorus(p.yzx, ab, r);
    vec4 d3 = sdgEllipticTorus(p.zxy, ab, r);

    d2.yzw = d2.wyz;
    d3.yzw = d3.zwy;

    vec4 res = (d1.x < d2.x) ? d1 : d2;
    res = (res.x < d3.x) ? res : d3;
    return res;
}

// 5. 0-Loop Branch-Cut-Free Conservative (p,q) Torus Knot SDF
float sdTorusKnot(vec3 p, vec2 torusR, vec2 pq, float rTube) {
    float R = torusR.x;
    float r1 = torusR.y;
    
    float rho = sqrt(p.x * p.x + p.z * p.z);
    float theta = atan(p.z, p.x);
    float u = sqrt((rho - R) * (rho - R) + p.y * p.y);
    float phi = atan(p.y, rho - R);
    
    float phase = pq.x * phi - pq.y * theta;
    float w = mod(phase + 3.141592653589793, 6.283185307179586) - 3.141592653589793;
    
    // First-order distance to the implicit helical phase locus
    float rhou = max(rho * u, 1e-10); // avoid div by zero
    float d_sheet = (rhou / sqrt(pq.x * pq.x * rho * rho + pq.y * pq.y * u * u)) * 2.0 * sin(abs(w) * 0.5);
    
    // Conservative chordal approximation
    float d_approx = sqrt(max(d_sheet * d_sheet + (u - r1) * (u - r1), 0.0)) - rTube;
    
    // Minor global factor for exactness across tight curvature regions
    float lipschitz = 0.85; 
    return d_approx * lipschitz;
}

#endif // SDF_GEMS_HAS_BORROMEAN_RINGS
