#ifndef SDF_GEMS_MILNOR_SEIFERT_KNOT_GLSL
#define SDF_GEMS_MILNOR_SEIFERT_KNOT_GLSL

// Complex power z^n = (r*e^{i*theta})^n for real exponent n >= 1.
vec2 sdfComplexPow(vec2 z, float n) {
    float r2 = dot(z, z);
    if (r2 < 1e-20) return vec2(0.0);
    float r = pow(r2, 0.5 * n);
    float ang = n * atan(z.y, z.x);
    return r * vec2(cos(ang), sin(ang));
}

// Pull back a 4D gradient G4 = (d/du_r, d/du_i, d/dv_r, d/dv_i) on S^3 subset C^2
// to 3D Euclidean space R^3 via the exact transpose Jacobian of inverse stereographic projection.
vec3 sdfPullbackStereoGrad(vec3 x, vec4 w, float mu, vec4 G4) {
    float dotWG = dot(w, G4);
    return mu * G4.xyz + (mu * (G4.w - dotWG)) * x;
}

// Evaluate Milnor polynomial f_{p,q}(u,v) = u^p - v^q on S^3 and its exact 3D spatial gradients
// gR = grad(Re f) and gI = grad(Im f).
vec2 sdfMilnorEval(vec3 p, float P, float Q, out vec3 gR, out vec3 gI) {
    float r2 = dot(p, p);
    float mu = 2.0 / (1.0 + r2);
    // w = (u_r, u_i, v_r, v_i) on S^3
    vec4 w = vec4(mu * p.x, mu * p.y, mu * p.z, 1.0 - mu);

    vec2 u = w.xy;
    vec2 v = w.zw;

    vec2 uP = sdfComplexPow(u, P);
    vec2 vQ = sdfComplexPow(v, Q);
    vec2 f = uP - vQ;

    // Complex Wirtinger derivatives: df/du = P * u^(P-1), df/dv = -Q * v^(Q-1)
    vec2 dU = P * sdfComplexPow(u, max(P - 1.0, 0.0));
    vec2 dV = -Q * sdfComplexPow(v, max(Q - 1.0, 0.0));

    // Cauchy-Riemann 4D gradients
    vec4 grad4R = vec4(dU.x, -dU.y, dV.x, -dV.y);
    vec4 grad4I = vec4(dU.y,  dU.x, dV.y,  dV.x);

    gR = sdfPullbackStereoGrad(p, w, mu, grad4R);
    gI = sdfPullbackStereoGrad(p, w, mu, grad4I);
    return f;
}

// Algebraic (P,Q)-Torus Knot SDF & analytical gradient via Milnor Fibration on S^3.
// Returns vec4(d_knot, grad_d_knot), plus writes solid-angle phase arg(f_{p,q}) in [-pi, pi] to outPhase.
vec4 sdgMilnorTorusKnot(vec3 p, float P, float Q, float tubeRadius, out float outPhase) {
    vec3 gR, gI;
    vec2 f = sdfMilnorEval(p, P, Q, gR, gI);
    outPhase = atan(f.y, f.x);

    float fMag = length(f);
    float gR2 = dot(gR, gR);
    float gI2 = dot(gI, gI);
    float gNorm = sqrt(0.5 * (gR2 + gI2) + 1e-10);

    // Distance estimate (scaled by 0.5 for conservative 1-Lipschitz bound across high (P,Q))
    float d = 0.5 * (fMag / gNorm) - tubeRadius;

    vec3 rawGrad = f.x * gR + f.y * gI;
    float rawLen = length(rawGrad);
    vec3 unitGrad = (rawLen > 1e-10) ? (rawGrad / rawLen) : normalize(gR + vec3(1e-6));
    return vec4(d, unitGrad);
}

// Branch-Cut-Free Seifert Minimal Soap-Film SDF bounded by the (P,Q)-Torus Knot at fiber angle alpha.
// Returns vec4(d_Seifert, grad_d_Seifert).
vec4 sdgSeifertSoapFilm(vec3 p, float P, float Q, float alpha, float thickness) {
    vec3 gR, gI;
    vec2 f = sdfMilnorEval(p, P, Q, gR, gI);

    float ca = cos(alpha), sa = sin(alpha);
    // Rotate f by e^{-i*alpha}
    float fR =  f.x * ca + f.y * sa;
    float fI = -f.x * sa + f.y * ca;
    vec3 tgR =  gR * ca + gI * sa;
    vec3 tgI = -gR * sa + gI * ca;

    float normR = sqrt(dot(tgR, tgR) + 1e-10);
    float normI = sqrt(dot(tgI, tgI) + 1e-10);

    // Normalized normal-plane coordinates (xi, eta) around the knot boundary
    float xi  = 0.5 * (fR / normR);
    float eta = 0.5 * (fI / normI);
    vec3 dirXi  = tgR / normR;
    vec3 dirEta = tgI / normI;

    // Distance to half-line {xi >= 0, eta == 0} in (xi, eta) plane (branch-cut-free!)
    float xiNeg = min(xi, 0.0);
    float distSheet = sqrt(xiNeg * xiNeg + eta * eta);
    float d = distSheet - thickness;

    vec3 grad = (xi >= 0.0)
        ? ((eta >= 0.0 ? 1.0 : -1.0) * dirEta)
        : normalize(xiNeg * dirXi + eta * dirEta + vec3(1e-10));
    return vec4(d, grad);
}

// Scalar convenience wrappers
float sdMilnorTorusKnot(vec3 p, float P, float Q, float tubeRadius) {
    float phase;
    return sdgMilnorTorusKnot(p, P, Q, tubeRadius, phase).x;
}

float sdSeifertSoapFilm(vec3 p, float P, float Q, float alpha, float thickness) {
    return sdgSeifertSoapFilm(p, P, Q, alpha, thickness).x;
}

#endif
