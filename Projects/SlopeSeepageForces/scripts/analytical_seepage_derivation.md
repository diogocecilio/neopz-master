# Optimal filtration velocity of the class (29): derivation and check of Eqs. 31-40

Reference: Ceron, Cecilio, Linn & Maghous, IJNAMG 2025, Sect. 4.1 (Eqs. 21-23, 29-40, Fig. 3).
Implementation: `analytical_seepage.py` (same folder). Run it to see all the checks quoted here.

## 0. Coordinates

Polar frame of Fig. 3: origin O at the crest edge, theta measured from the crest (negative x axis)
towards the face. In the shared paper coordinates (x to the right, y DOWN):

    x = -r cos(theta),  y = r sin(theta),  theta = atan2(y, -x)
    e_r = (-cos theta, sin theta),   e_theta = (sin theta, cos theta)
    e_y = sin(theta) e_r + cos(theta) e_theta

- crest: theta = 0, outward normal n = -e_y = -e_theta, u^d = 0
- face: theta = Theta := pi - beta, outward normal n = e_theta = (sin beta, -cos beta),
  u^d = -gamma_w r sin(beta) for r < R_w (above the lowered water), u^d = -gamma_w h_w for R_w < r < R
- toe ground: y = H, i.e. theta_m(r) = pi - arcsin(H/r) for r > R, n = -e_y, u^d = -gamma_w h_w
- R_w = h_w/sin(beta), R = H/sin(beta), R_e = sqrt(H^2 + (L_m + H/tan beta)^2)

K^-1 = (1/k_h) [ 1 + (alpha - 1) e_y (x) e_y ], so with v = v_r e_r + v_t e_theta

    k_h v.K^-1.v = c(theta) v_r^2 + d(theta) v_t^2 + (alpha-1) sin(2 theta) v_r v_t
    c = cos^2 + alpha sin^2,   d = sin^2 + alpha cos^2.

J*(v') = 1/2 int v'.K^-1.v' dOmega + int_{dOmega_u} u^d v'.n dS (Eq. 23). The crest carries no
boundary term (u^d = 0). Fields of the class satisfy v.e_r = 0 at r = R_w, R, R_e, so the zones do
not exchange flux and J* = J1*(h1, h2) + J2*(h3) + J3*(h4) (Eq. 30): each zone is optimised separately.

## 1. Zones 2 and 3 (Eqs. 32-35): CONFIRMED

Zone 2, v = h3(r) e_theta: int_0^Theta d dtheta = (alpha+1)(pi-beta)/2 - (alpha-1) sin(2 beta)/4 = A(alpha)
(Eq. 34), face term -gamma_w h_w int h3 dr, so

    J2* = int_{R_w}^{R} [ A r h3^2 / (2 k_h) - gamma_w h_w h3 ] dr
    => h3 = k_h gamma_w h_w / (A r)          (Eq. 32),   J2* = -(k_h gamma_w^2 h_w^2 / (2A)) ln(H/h_w).

Zone 3, v = h4(r) e_theta on 0 < theta < theta_m(r): int_0^{theta_m} d dtheta = B(r, alpha)/2 with B of
Eq. 35 (sin 2 theta_m = -2 H sqrt(r^2-H^2)/r^2). On the toe ground v.n = -h4 cos(theta_m) =
h4 sqrt(r^2 - H^2)/r and dS = dX = r dr / sqrt(r^2 - H^2), so v.n dS = h4 dr exactly. Hence

    J3* = int_R^{R_e} [ B r h4^2 / (4 k_h) - gamma_w h_w h4 ] dr
    => h4 = 2 k_h gamma_w h_w / (B r)        (Eq. 33),   J3* = -(k_h gamma_w^2 h_w^2 / 4) int_R^{R_e} 4 dr/(B r).

## 2. Zone 1: optimum in h1 for a given h2

Put g(r) = r h1(r) (stream function psi = -g(r) h2(theta)), so v_r = -(g/r) h2', v_t = g' h2, g(R_w) = 0
and g(0) = 0 (finite energy).

- The cross term integrates to (alpha-1) int sin2theta h2 h2' dtheta * int g g' dr = ... * [g^2/2]_0^{R_w} = 0.
- Face term: int_0^{R_w} (-gamma_w r sin beta) g' h2(Theta) dr = + gamma_w sin(beta) h2e int_0^{R_w} g dr
  (integration by parts, g(0) = g(R_w) = 0), with h2e := h2(Theta).

With C, D of Eq. 36 (C = int c h2'^2, D = int d h2^2):

    J1* = (1/(2 k_h)) [ C int g^2/r dr + D int r g'^2 dr ] + gamma_w sin(beta) h2e int g dr.

Euler-Lagrange in g:   D (r g')' - C g / r = k_h gamma_w sin(beta) h2e.
With m := sqrt(C/D), the solution regular at 0 and vanishing at R_w is (s = r/R_w)

    g = a r [1 - s^(m-1)],   a = k_h gamma_w sin(beta) h2e / (D - C),
    v1 = [k_h gamma_w sin(beta) h2e / (D - C)] { [1 - m s^(m-1)] h2 e_theta - [1 - s^(m-1)] h2' e_r }.   (31-corrected)

**Disagreement with the printed Eq. 31:**
1. prefactor: printed 1/(C - D); the correct sign is 1/(D - C). With the printed sign the face flux
   above the water is reversed and J* > 0 (check: beta=30, alpha=1, h_w=H gives J* = +48.9 instead of
   J*(opt) = -41.6, k_h = 1, gamma_w = 9.81, H = 1).
2. e_theta term: printed (r/R_w)^(C/D - 1) * C/D; it must be (r/R_w)^(sqrt(C/D) - 1) * sqrt(C/D), the same
   exponent as the e_r term, because v_t must equal (r h1)' for div v = 0. The printed field is not
   divergence-free (relative residual 0.64 in that example) and is therefore not admissible.
Both are typesetting errors: the rest of the paper (Eqs. 37, 40 and Fig. 5) is consistent with the corrected form.

Since a(1 - m) = k_h gamma_w sin(beta) h2e / (D(1 + m)) =: a', the code uses the form that stays regular at m = 1:
v_r = -a' E h2', v_t = a' (E + s^(m-1)) h2, with E = (1 - s^(m-1))/(1 - m) (E -> ln s as m -> 1).

Value at the optimum g. With G0 = int g^2/r = a^2 R_w^2 (m-1)^2/(2m(m+1)), G1 = int r g'^2 =
a^2 R_w^2 (m-1)^2/(2(m+1)), G2 = int g = a R_w^2 (m-1)/(2(m+1)) and J1* = (1/2) gamma_w sin(beta) h2e G2
(quadratic functional at its minimum):

    J1* = -(k_h gamma_w^2 h_w^2 / 4) h2e^2 / (sqrt C + sqrt D)^2        (first term of Eq. 40: CONFIRMED)

## 3. Optimum in h2 (Eqs. 36-39)

J1* depends on h2 only through F[h2] = h2e^2 / (sqrt C + sqrt D)^2, which must be MAXIMISED (homogeneous of
degree 0: (h1, h2) -> (lambda h1, h2/lambda) leaves v unchanged).

Stationarity (equivalently, the first variation of J1*(g, h2) in h2 at the optimal g, where G1/G0 = m):

    (c h2')' - m d h2 = 0,  i.e.  h2'' + (alpha-1) sin2theta / c * h2' - sqrt(C/D) * d/c * h2 = 0   (Eq. 37: CONFIRMED)
    h2'(0) = 0                                                                                    (Eq. 38: CONFIRMED)
    c(Theta) h2'(Theta) h2(Theta) = sqrt(C) (sqrt C + sqrt D)                                     (natural BC at the face)

The sqrt(C/D) coefficient of Eq. 37 is correct (it is the ratio G1/G0 = m of the optimal g), not a typo.

Multiplying (37) by h2 and integrating shows that the face condition above is an IDENTITY for every
solution of (37)-(38) with m = sqrt(C/D) (checked numerically to 1e-13). The printed Eq. 39, after
(D - C)^2/(sqrt C - sqrt D)^2 = (sqrt C + sqrt D)^2, reads
h2'(Theta) = 2 sqrt C (sqrt C + sqrt D)^3 / (R_w^2 c(Theta) k_h gamma_w sin(beta) h2e^2), i.e. the identity times
the factor 2 (sqrt C + sqrt D)^2 / (R_w^2 k_h gamma_w sin(beta) h2e) = 1. It is not homogeneous in h2, so it
only fixes the free amplitude lambda of h2 (h2e = 2 (sqrt C + sqrt D)^2 / (R_w^2 k_h gamma_w sin beta)).
**Printed Eq. 39 is consistent but is a normalisation, not an extra condition; v'_opt, f and J* do not depend on it.**
The code normalises h2(0) = 1 instead.

Robust solution (no fixed-point iteration on C, D needed). Cauchy-Schwarz gives
(sqrt C + sqrt D)^2 = min_{t>0} (1 + t)(C + D/t) (equality at t = 1/m), hence

    1/F_max = min_m Phi(m),  Phi(m) = (1 + 1/m) min_{h2(Theta)=1} (C + m D) = (1 + 1/m) c(Theta) phi'(Theta)/phi(Theta),

with phi the solution of the LINEAR problem (c phi')' = m d phi, phi(0) = 1, phi'(0) = 0 (LSODA). The
stationarity of Phi in m is exactly m = sqrt(C/D) (envelope theorem). The code scans ln m, refines with
Brent, sharpens the root of m - sqrt(C(m)/D(m)) (residual 1e-11), and cross-checks against a direct
L-BFGS minimisation of (sqrt C + sqrt D)^2 over a 400-element P1 discretisation of h2 (no ODE),
which agrees to about 1e-5 relative.

Degenerate case. Phi(0+) = A. For alpha = 1, phi = cosh(sqrt(m) theta) and
Phi(m) = Theta + m Theta (1 - Theta^2/3) + O(m^2): for Theta < sqrt 3 (beta > 80.76 deg) Phi increases from
m = 0 and the optimum is the limit m -> 0 of the class: h2 = const, purely tangential zone-1 field
v = (k_h gamma_w sin(beta)/A) e_theta, F = 1/A. Numerical thresholds: alpha = 1: 80.76 deg; alpha = 2: 81.10 deg;
alpha = 4: 88.40 deg; alpha >= 5: none up to 90 deg. The paper's Fig. 5 values agree with this limit
(e.g. alpha = 1, beta = 90: 0.5914 in both), so the authors also reached it.

## 4. Final functional (Eq. 40): CONFIRMED

    J*(v'_opt) = -(k_h h_w^2 gamma_w^2 / 4) [ h2e^2/(sqrt C + sqrt D)^2 - (2/A) ln(h_w/H) + int_R^{R_e} 4 dr/(B r) ].

Fig. 5 (h_w = H, L_m = 10 H): the solid curves extracted from the vector PDF (`data/fig5_vector_fill_polygons.csv`)
are reproduced at all 64 points (alpha = 1, 2, 4, 10; beta = 15..90 by 5 deg) within 0.02 %
(`results/analytical_seepage/fig5_analytical_Jstar.csv`).

## 5. Remarks on the resulting field (relevant for the stability step)

- No flux crosses r = R_w, R, R_e. Inside r < R_e, v'_opt does not depend on L_m; only J* does, through
  zone 3: J* ~ -(k_h gamma_w^2 h_w^2/((alpha+1) pi)) ln R_e (log-divergent, as the paper notes).
- Zone 1 is a closed flow cell: psi = -g h2 vanishes at r = 0 and r = R_w, so the net flux through the
  upper face (and through the crest segment 0 < -x < R_w) is zero. For m != 1 (either side of 1), v.n < 0 (inflow) on the
  face for r < s* R_w with s* = m^(1/(1-m)) (0.34 for beta = 30, alpha = 1), and the matching outflow leaves
  through the crest near O: there is a small recirculation cell at the crest edge where f = K^-1 v points
  up and away from the face. Elsewhere f points down, turns towards the face and exits upward
  through the toe ground, as in Fig. 4c.
- For 0 < m < 1, |v| ~ r^(m-1) near O (integrable singularity). The exact solution is regular there
  (locally u = -gamma_w y, f = gamma_w e_y).
- The tangential component of v jumps across r = R_w and r = R. This is allowed, because only v.n has to be continuous.
- h_w = 0 gives v = 0 (`AnalyticalSeepage(..., hw=0).force == 0`).
- Degenerate case: the m -> 0 field (h2 = const, h1 = a (1 - R_w/r)) is admissible in V (div-free, bounded,
  normal flux continuous) but h1 is unbounded at r = 0, so it is the infimum of J* over the class (29), not a
  member of it. Fields of the class with g(0) = 0 approach it from above only logarithmically (see check C).

## 6. Independent (adversarial) check: `check_analytical_seepage.py`

Code written separately from `analytical_seepage.py` (only the module's public API is called).
Log: `results/analytical_seepage/independent_check_output.txt`.

- A) The h2 problem was solved again by a Chebyshev-Ritz method (no ODE). F matches the module to 1e-13. m matches to 1e-8.
  The degenerate/non-degenerate flag is the same for every beta from 15 to 90 deg (1 deg steps) at alpha = 1, 5 and 10 (option `--sweep`).
  The small-m expansion Phi(m) = A + m (A - P) + O(m^2), with P = int_0^Theta (int_0^t d)^2 / c dt, gives the
  threshold angles A = P: 80.761 deg (alpha = 1; analytically pi - sqrt 3), 81.101 deg (alpha = 2), 88.398 deg (alpha = 4), and none below 90 deg for alpha >= 5.
- B) J* of the implemented field was recomputed by a 2-D quadrature of Eq. 23 (Cartesian v, Cartesian K^-1,
  Cartesian normals, and the toe-ground angle from atan2 instead of Eq. 35). It agrees with Jstar() to 1e-12 or better.
- C) J* was minimised directly with BFGS over a discretised class:
  - g = r h1 as a graded C1 cubic Hermite function, with g(0) = g(R_w) = 0;
  - h2, h3 as Chebyshev series, and r h4 as a Chebyshev series in ln(x + H).

  Results in the non-degenerate cases: J*_direct lies 2e-6 to 4e-6 (relative) above Jstar() and decreases further as the mesh is refined.
  The velocity matches to 1e-8 in zones 2 and 3, and to (1-3)e-3 of max|v| in zone 1 (discretisation error, which shrinks with refinement).

  Results in the degenerate case (beta = 85, alpha = 1): with h2 = const and g(0) free, J* is reproduced to 1e-15. Inside the class
  (g(0) = 0), J*_direct - J*(m -> 0) is about 0.072 |J*| / ln(1/s_min). It is positive and tends to zero: no discretised member of the class beats the m -> 0 value, as the analytical reduction predicts.
- D) The printed Eq. 31 gives J* = +48.88 against -41.62 at the optimum (beta = 30, alpha = 1, h_w = H), and div v != 0.
  Fixing only the sign gives J* = -48.03, which is BELOW the optimum. That is possible only because the field is not
  divergence-free, so it is not admissible. Both typos must therefore be corrected.
