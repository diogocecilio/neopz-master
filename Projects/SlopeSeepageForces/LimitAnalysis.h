// Kinematic limit analysis of Ceron et al. (IJNAMG 2025), Sect. 4.2 (Figs. 6-7, Eqs. 42-58): upper bound of the
// stability factor Gamma = min P_mr / (P_gamma + P_u) (Eq. 55) over rotational log-spiral mechanisms, with the
// seepage forces of any slope::ForceField. C++ port of scripts/limit_analysis.py (the verified reference): same
// parameterization, exact admissibility test, closed forms, quadrature rules and optimiser settings.
//
// Frame: internally the paper's frame of the Python reference (origin O at the crest edge, x to the right, y DOWN,
// toe T = (H / tan(beta), H)), so that every formula reads as in the paper and in limit_analysis.py. The seepage
// fields come in NeoPZ coordinates (y up) and are converted at each call: point (x, y_paper) -> (x, -y_paper),
// force (fx, fy) -> (fx, -fy), u unchanged (f . U is invariant under the reflection). Results give A, B, C in both.
//
// Mechanisms (angles measured at the centre C from the -x direction, turning downwards, Fig. 7):
//   point(rho, theta) = C + rho e_r, e_r = (-cos theta, sin theta), e_theta = (sin theta, cos theta), C = (Cx, -Cy)
//   (Cy = height of C above the crest); log spiral r(theta) = r0 exp((theta - theta1) tan phi) from A on the crest
//   (theta1) to B on the ground surface (theta2); rigid rotation U = omega (y + Cy, Cx - x), omega = 1 below.
//   Mechanism I : B on the face, B = (eta H / tan beta, eta H), 0 < eta <= 1 (Fig. 7);
//   Mechanism II: B on the toe ground, B = (H / tan beta + d, H), 0 <= d <= d_max (default 10 H, Fig. 6 right).
//   Unified parameter s = eta (I) or 1 + d / H (II). Given (theta1, theta2, s):
//     r0 = B_y / (E sin theta2 - sin theta1), E = exp((theta2 - theta1) tan phi), Cx = B_x + r0 E cos theta2,
//     Cy = r0 sin theta1, L = OA = r0 cos theta1 - Cx (Eqs. 43-46).
//   Admissibility (Eq. 56 and the geometric conditions it implies): 0 < theta1 < theta2 < pi, r0 > 0 with a
//   denominator not dominated by round-off, r_h <= 1000 H, L > 0 (A on the crest left of O), C on the air side of the
//   face line, theta_O in [theta1, theta_F], theta2 < pi - beta (I), spiral below O (and below T with
//   theta_T <= theta2 for II). Along each straight piece of the surface ln(r(theta) / rho_surface(theta)) is concave,
//   so these point tests are exact.
// Rates of work per unit omega:
//   P_mr    = c r0^2 (exp(2 (theta2 - theta1) tan phi) - 1) / (2 tan phi)  (Eq. 48; c r0^2 (theta2 - theta1), phi = 0)
//   P_gamma = gamma' (f1 - f2 - f3 [- f4 for II]) (Eqs. 50-53, f3 in a form valid at beta = 90 deg, f4 = triangle
//             C-T-B of mechanism II), gamma' = gamma - gamma_w; also by quadrature (both mechanisms);
//   P_u     = int_Omega f . U dOmega (Eq. 49), either
//             (domain)   composite Gauss rule in polar coordinates about C: per ray theta, rho from the ground surface
//                        to the spiral; angular range split at the rays through O and T, panels graded towards the
//                        piece ends; rays split where they cross the circles of the field, and angular breaks at the
//                        rays tangent to its discontinuity circles and through their intersections with the surface and
//                        the spiral (analytical field: jumps on r = R_w, R, R_e about O). Jumps on other curves (FE
//                        element edges) are not split: the error then decreases algebraically with the level;
//             (boundary) for f = -grad u (FE fields; exact for the discrete u_h): P_u = -closed int u U.n ds since
//                        div U = 0 for a rigid rotation: tan(phi) int u r^2 dtheta on the spiral (U.n = -|U| sin phi)
//                        minus the ground-surface terms on A-O, O-F and T-B (u only on the boundary).
// Optimisation: for each mechanism class and seed, global-best PSO (constriction coefficients of Clerc & Kennedy,
// 40 particles, <= 150 iterations, stall 40) started from admissible random points with the coarse quadrature (more
// random points than the Python when mechanisms with P_ext > 0 are rare, see PSO), then a
// bounded Nelder-Mead polish (scipy's algorithm, 2 restarts, xatol 1e-10, fatol 1e-13) with the fine quadrature;
// Gamma = best over the runs. The random streams are those of std::mt19937_64 (not numpy's), so the PSO paths differ
// from the Python ones while the polished optima agree to the optimiser tolerance. The runs (class x seed) are
// independent and are spread over threads; the result does not depend on the number of threads (fields are read-only).
// Gamma = inf (BIG) when no admissible mechanism has P_gamma + P_u > 0 (e.g. beta < phi without seepage).
// By similarity H_crit = Gamma(H) H at fixed h_w / H, beta, alpha, phi (see SPEC): P_mr ~ c H^2, P_gamma, P_u ~ gamma H^3.
#ifndef LIMITANALYSIS_H
#define LIMITANALYSIS_H

#include "SeepageForceField.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <thread>
#include <vector>

namespace slope {
namespace la {

constexpr REAL BIG = 1.e30;     ///< objective value of inadmissible mechanisms (finite: clean Nelder-Mead arithmetic)
constexpr REAL EPS_ANG = 1.e-6; ///< margin of the angular search bounds
constexpr REAL R_MAX = 1000.;   ///< admissible mechanisms: r_h <= R_MAX H (beyond, the closed forms lose all digits)

/// excess pore pressure u at x (NeoPZ coordinates); false outside the domain of the field
using PotentialField = std::function<bool(const TPZVec<REAL> &x, REAL &u)>;

/// u of a FE field (the P2 interpolant of PoreField; f = -grad u exactly)
inline PotentialField PotentialOf(std::shared_ptr<const PoreField> pf) {
    return [pf](const TPZVec<REAL> &x, REAL &u) {
        REAL g[2];
        return pf->Evaluate(x, u, g);
    };
}

/// circle (NeoPZ coordinates) where the field may jump or is singular: the rays of the domain quadrature are split
/// where they cross it; a discontinuity circle also adds angular breaks (see DomainPower)
struct Circle {
    REAL x = 0., y = 0., R = 0.;
    bool discontinuity = true;
};

/// seepage input of the limit analysis (NeoPZ coordinates)
struct Seepage {
    ForceField force;             ///< f (kN/m^3); empty with u empty: no seepage (h_w = 0, P_u = 0)
    PotentialField u;             ///< u with f = -grad u (optional): enables the boundary formula
    bool gradient = false;        ///< f = -grad u exactly (FE fields): method "auto" uses the boundary formula
    REAL hw = -1.;                ///< water line depth below the crest (m): panel break on the face (boundary formula)
    std::vector<Circle> circles;  ///< curves where f jumps (analytical field: r = R_w, R, R_e about O)
    Circle grade{0., 0., 0., false}; ///< grade.R > 0: circles grade.R gradeRatio^k, k = 1 .. Quadrature::nCircleGrade,
    REAL gradeRatio = 0.3;           ///< about (grade.x, grade.y): grading towards a singular point (r^(m-1) at O)
    bool None() const { return !force && !u; }
};

/// FE field f = -grad u_h of the seepage solution, water line at depth hw (m)
inline Seepage FESeepage(std::shared_ptr<const PoreField> pf, REAL hw) {
    Seepage s;
    s.force = pf->AsForceField(pf);
    s.u = PotentialOf(pf);
    s.gradient = true;
    s.hw = hw;
    return s;
}

/// slope geometry and friction angle (the soil strength enters only through c in P_mr)
struct Geometry {
    REAL H = 1., betaDeg = 45., phiDeg = 30.;
    REAL beta = 0., sb = 0., cb = 0., xtoe = 0., k = 0.; ///< beta (rad), sin, cos (0 at 90 deg), x of T, tan(phi)
    Geometry() = default;
    Geometry(REAL H_, REAL betaDeg_, REAL phiDeg_) : H(H_), betaDeg(betaDeg_), phiDeg(phiDeg_) {
        beta = betaDeg * M_PI / 180.;
        sb = std::sin(beta);
        cb = std::fabs(betaDeg - 90.) < 1.e-12 ? 0. : std::cos(beta);
        xtoe = H * cb / sb;
        k = std::tan(phiDeg * M_PI / 180.);
    }
};

/// one log-spiral mechanism (paper coordinates, see the header)
struct Mechanism {
    REAL t1 = 0., t2 = 0., s = 0.; ///< theta1, theta2, s
    bool II = false;               ///< B on the toe ground (s > 1)
    bool ok = false;               ///< admissible
    REAL Bx = 0., By = 0., E = 0., den = 0., r0 = 0., rh = 0., Cx = 0., Cy = 0., Ax = 0., L = 0.;
    REAL thO = 0., dO = 0., thT = 0., dT = 0., D = 0., thF = 0.; ///< rays through O and T, distance C - face line, end of the face piece
    REAL R(REAL theta, REAL k) const { return r0 * std::exp((theta - t1) * k); }
};

inline Mechanism MakeMechanism(const Geometry &g, REAL t1, REAL t2, REAL s) {
    Mechanism m;
    m.t1 = t1, m.t2 = t2, m.s = s;
    const REAL k = g.k, H = g.H, xt = g.xtoe;
    m.II = s > 1.;
    m.Bx = m.II ? xt + (s - 1.) * H : s * xt;
    m.By = m.II ? H : s * H;
    m.E = std::exp((t2 - t1) * k);
    m.den = m.E * std::sin(t2) - std::sin(t1);
    m.r0 = m.By / m.den;
    m.rh = m.r0 * m.E;
    m.Cx = m.Bx + m.rh * std::cos(t2);
    m.Cy = m.r0 * std::sin(t1);
    m.Ax = m.Cx - m.r0 * std::cos(t1);
    m.L = -m.Ax;
    m.thO = std::atan2(m.Cy, m.Cx); // Eq. 47
    m.dO = std::hypot(m.Cx, m.Cy);
    m.thT = std::atan2(m.Cy + H, m.Cx - xt);
    m.dT = std::hypot(m.Cx - xt, m.Cy + H);
    m.D = m.Cx * g.sb + m.Cy * g.cb;
    m.thF = m.II ? m.thT : t2;
    bool ok = t1 > 0. && t2 > t1 && t2 < M_PI && s > 0.;
    ok = ok && m.den > 1.e-12 * std::max(m.E * std::fabs(std::sin(t2)), 1.e-300); // no round-off r0
    ok = ok && std::isfinite(m.r0) && m.r0 > 0. && m.rh <= R_MAX * H;
    ok = ok && m.L > 0. && m.D > 0.;
    ok = ok && m.thO >= t1 && m.thO <= m.thF;
    if (!m.II) ok = ok && t2 < M_PI - g.beta + 1.e-12; // Eq. 56
    if (ok) {
        const REAL lr0 = std::log(m.r0);
        ok = (m.thO - t1) * k + lr0 >= std::log(m.dO) - 1.e-12;                                  // spiral below O
        if (ok && m.II) ok = (m.thT - t1) * k + lr0 >= std::log(m.dT) - 1.e-12 && m.thT <= t2; // spiral below T
    }
    m.ok = ok;
    return m;
}

/// Eq. 48 divided by omega
inline REAL Pmr(const Mechanism &m, const Geometry &g, REAL c) {
    const REAL dth = m.t2 - m.t1, k = g.k;
    const REAL fac = k < 1.e-14 ? dth : std::expm1(2. * dth * k) / (2. * k);
    return c * m.r0 * m.r0 * fac;
}

/// f1, f2, f3 (Eqs. 51-53, f3 in a form valid for beta = 90 deg) and f4 (triangle C-T-B, mechanism II)
inline void FTerms(const Mechanism &m, const Geometry &g, REAL f[4]) {
    const REAL k = g.k, r0 = m.r0, rh = m.rh, t1 = m.t1, t2 = m.t2;
    f[0] = (rh * rh * rh * (3. * k * std::cos(t2) + std::sin(t2)) - r0 * r0 * r0 * (3. * k * std::cos(t1) + std::sin(t1))) /
           (3. * (9. * k * k + 1.));
    const REAL s1 = std::sin(t1), sO = std::sin(m.thO);
    f[1] = r0 * r0 * r0 * s1 / 6. * (1. - s1 * s1 / (sO * sO));
    auto G = [&g](REAL psi) {
        const REAL sp = std::sin(psi);
        return -g.cb / (2. * sp * sp) - g.sb * std::cos(psi) / sp;
    };
    f[2] = m.D * m.D * m.D / 3. * (G(m.thF + g.beta) - G(m.thO + g.beta));
    if (m.II) {
        const REAL a = g.H + m.Cy, sT = std::sin(m.thT), s2 = std::sin(t2);
        f[3] = a * a * a / 6. * (1. / (sT * sT) - 1. / (s2 * s2));
    } else
        f[3] = 0.;
}

/// Eq. 53 as printed (mechanism I, beta < 90 deg), for the checks
inline REAL F3Paper(const Mechanism &m, const Geometry &g) {
    const REAL tb = std::tan(g.beta), a = m.Cy + tb * m.Cx;
    const REAL d1 = std::tan(m.thO) + tb, d2 = std::tan(m.t2) + tb;
    return a * a * a / 6. * (1. / (d1 * d1) - 1. / (d2 * d2));
}

/// Eq. 50 (+ f4 for mechanism II) divided by omega
inline REAL PgammaClosed(const Mechanism &m, const Geometry &g, REAL gammap) {
    REAL f[4];
    FTerms(m, g, f);
    return gammap * (f[0] - f[1] - f[2] - f[3]);
}

/// Gauss-Legendre rule with q points on [0, 1] (ascending nodes)
inline void Gauss01(int q, std::vector<REAL> &x, std::vector<REAL> &w) {
    x.assign(q, 0.), w.assign(q, 0.);
    for (int i = 0; i < q; i++) {
        REAL z = std::cos(M_PI * (i + 0.75) / (q + 0.5)), dp = 1.;
        for (int it = 0; it < 100; it++) {
            REAL p0 = 1., p1 = z; // Legendre recurrence
            for (int n = 2; n <= q; n++) {
                const REAL p2 = ((2. * n - 1.) * z * p1 - (n - 1.) * p0) / n;
                p0 = p1, p1 = p2;
            }
            const REAL pq = q == 1 ? z : p1, pm = q == 1 ? 1. : p0;
            dp = q * (z * pq - pm) / (z * z - 1.);
            const REAL dz = pq / dp;
            z -= dz;
            if (std::fabs(dz) < 1.e-16) break;
        }
        // node z (descending in i): store ascending on [0, 1]
        x[q - 1 - i] = 0.5 * (1. + z);
        w[q - 1 - i] = 1. / ((1. - z * z) * dp * dp); // = 0.5 * 2 / ((1 - z^2) P'^2)
    }
}

/// Composite Gauss rule in polar coordinates about C (limit_analysis.py Quadrature):
///   theta: nPanels uniform panels per piece + nGrade geometric panels (ratio gradeRatio) at both ends of each piece
///          (rays through A, O, T, B), qTheta points per panel;
///   rho  : nRho uniform panels between the ground surface and the spiral, further split at the intersections with
///          the field circles, qRho points per sub-segment; nCircleGrade = number of grading circles (Seepage::grade).
struct Quadrature {
    std::string name;
    int nPanels = 6, qTheta = 6, nGrade = 2, nRho = 2, qRho = 5, nCircleGrade = 0;
    REAL gradeRatio = 0.25;
    std::vector<REAL> br;     ///< panel breaks of one piece on [0, 1]
    std::vector<REAL> xt, wt; ///< angular Gauss rule on [0, 1]
    std::vector<REAL> xr, wr; ///< radial Gauss rule on [0, 1]

    Quadrature() { Init(); }
    Quadrature(const std::string &nm, int np, int qt, int ng, int nr, int qr, int ncg, REAL ratio = 0.25)
        : name(nm), nPanels(np), qTheta(qt), nGrade(ng), nRho(nr), qRho(qr), nCircleGrade(ncg), gradeRatio(ratio) {
        Init();
    }
    void Init() {
        const REAL h = 1. / nPanels;
        std::set<REAL> b;
        for (int i = 0; i <= nPanels; i++) b.insert(i == nPanels ? 1. : i * h);
        for (int j = 1; j <= nGrade; j++) {
            b.insert(h * std::pow(gradeRatio, j));
            b.insert(1. - h * std::pow(gradeRatio, j));
        }
        br.assign(b.begin(), b.end());
        Gauss01(qTheta, xt, wt);
        Gauss01(qRho, xr, wr);
    }
    /// levels of limit_analysis.py (QUAD_LEVELS): coarse (PSO search), medium, fine (final), xfine, ref
    static Quadrature Level(const std::string &nm) {
        if (nm == "coarse") return Quadrature(nm, 4, 5, 1, 1, 4, 0);
        if (nm == "medium") return Quadrature(nm, 8, 6, 3, 2, 5, 2);
        if (nm == "fine") return Quadrature(nm, 12, 8, 5, 2, 6, 4);
        if (nm == "xfine") return Quadrature(nm, 24, 10, 8, 4, 8, 7);
        if (nm == "ref") return Quadrature(nm, 48, 12, 12, 6, 10, 10);
        if (nm == "dense") return Quadrature(nm, 120, 6, 6, 60, 4, 10); // check 8 of limit_analysis.py (FE fields)
        return Quadrature("", 0, 0, 0, 0, 0, 0);                           // nPanels = 0: unknown name
    }
    bool Valid() const { return nPanels > 0; }
    /// panels and points of the boundary formula: 24 x 8 for coarse / medium, 48 x 8 otherwise (as limit_analysis.py)
    int BoundaryPanels() const { return (name == "coarse" || name == "medium") ? 24 : 48; }
};

/// circle in paper coordinates (qx, qy_paper, R)
struct PCircle {
    REAL qx, qy, R;
    bool disc;
};

/// Angles (seen from C) where the angular integrand has kinks because of the circle: rays tangent to the circle,
/// rays through its intersections with the ground surface (<= 6) and with the spiral (<= 4, sign changes of
/// |S(theta) - Q|^2 - R^2 on 48 samples refined by 48 bisections). out[12], NaN where absent.
inline void CircleAngles(const Mechanism &m, const Geometry &g, const PCircle &c, REAL out[12]) {
    const REAL nan = std::numeric_limits<REAL>::quiet_NaN();
    const REAL qx = c.qx, qy = c.qy, R = c.R, H = g.H, xt = g.xtoe, sb = g.sb, cb = g.cb;
    const REAL Cx = m.Cx, Cy = m.Cy;
    for (int i = 0; i < 12; i++) out[i] = nan;
    // tangent rays
    const REAL dx = Cx - qx, dy = Cy + qy, dist = std::hypot(dx, dy), thc = std::atan2(dy, dx);
    if (dist > R) {
        const REAL a = std::asin(R / dist);
        out[0] = thc - a, out[1] = thc + a;
    }
    // intersections with the ground surface (fixed points)
    REAL pts[6][2];
    int np = 0;
    if (R > std::fabs(qy)) {
        const REAL h = std::sqrt(R * R - qy * qy);
        for (REAL x0 : {qx - h, qx + h})
            if (x0 <= 0. && np < 6) pts[np][0] = x0, pts[np][1] = 0., np++;
    }
    const REAL dq = qx * cb + qy * sb, disc = dq * dq - (qx * qx + qy * qy - R * R);
    if (disc > 0.) {
        const REAL sq = std::sqrt(disc);
        for (REAL t : {dq - sq, dq + sq})
            if (t >= 0. && t <= H / sb && np < 6) pts[np][0] = t * cb, pts[np][1] = t * sb, np++;
    }
    if (R > std::fabs(H - qy)) {
        const REAL h = std::sqrt(R * R - (H - qy) * (H - qy));
        for (REAL x0 : {qx - h, qx + h})
            if (x0 >= xt && np < 6) pts[np][0] = x0, pts[np][1] = H, np++;
    }
    for (int i = 0; i < np; i++) out[2 + i] = std::atan2(pts[i][1] + Cy, Cx - pts[i][0]);
    // intersections with the spiral
    const REAL k = g.k, t1 = m.t1, t2 = m.t2, r0 = m.r0;
    auto gfun = [&](REAL th) {
        const REAL r = r0 * std::exp((th - t1) * k);
        const REAL ex = Cx - r * std::cos(th) - qx, ey = -Cy + r * std::sin(th) - qy;
        return ex * ex + ey * ey - R * R;
    };
    constexpr int ns = 48, nb = 48;
    REAL ts[ns], gs[ns];
    for (int i = 0; i < ns; i++) {
        ts[i] = t1 + (t2 - t1) * (i == ns - 1 ? 1. : i * (1. / (ns - 1)));
        gs[i] = gfun(ts[i]);
    }
    int nr = 0;
    for (int j = 0; j < ns - 1 && nr < 4; j++) {
        if (std::signbit(gs[j]) == std::signbit(gs[j + 1])) continue;
        REAL lo = ts[j], hi = ts[j + 1], glo = gs[j];
        for (int it = 0; it < nb; it++) {
            const REAL mid = 0.5 * (lo + hi), gm = gfun(mid);
            if (std::signbit(gm) == std::signbit(glo)) lo = mid, glo = gm;
            else hi = mid;
        }
        out[8 + nr++] = 0.5 * (lo + hi);
    }
}

/// circles of the field for a quadrature level, in paper coordinates
inline std::vector<PCircle> FieldCircles(const Seepage &sp, int nCircleGrade) {
    std::vector<PCircle> out;
    for (const Circle &c : sp.circles)
        if (c.R > 0.) out.push_back({c.x, -c.y, c.R, c.discontinuity});
    if (sp.grade.R > 0.)
        for (int k = 1; k <= nCircleGrade; k++) out.push_back({sp.grade.x, -sp.grade.y, sp.grade.R * std::pow(sp.gradeRatio, k), false});
    return out;
}

/// int_Omega f . U dOmega / omega for an admissible mechanism; fpaper(x, y_paper, f_paper[2]) (see the header for
/// the rule). The integrand is rho^2 (f_x sin theta + f_y cos theta) (U = rho e_theta, dOmega = rho drho dtheta).
template <class F>
REAL DomainPower(const Mechanism &m, const Geometry &g, const Quadrature &q, const F &fpaper,
                 const std::vector<PCircle> &circles) {
    const REAL t1 = m.t1, t2 = m.t2, k = g.k, H = g.H;
    const REAL ta[3] = {t1, m.thO, m.II ? m.thT : t2};
    const REAL tb[3] = {std::max(m.thO, ta[0]), std::max(m.thF, ta[1]), std::max(t2, ta[2])};
    std::vector<REAL> brk;
    brk.reserve(3 * q.br.size() + 12 * circles.size());
    for (int p = 0; p < 3; p++)
        for (REAL b : q.br) brk.push_back(ta[p] + (tb[p] - ta[p]) * b);
    for (const PCircle &c : circles) {
        if (!c.disc) continue;
        REAL sa[12];
        CircleAngles(m, g, c, sa);
        for (REAL a : sa) brk.push_back(std::isfinite(a) ? std::min(std::max(a, t1), t2) : t1);
    }
    std::sort(brk.begin(), brk.end());
    std::vector<REAL> rb;
    rb.reserve(q.nRho + 1 + 2 * circles.size());
    REAL total = 0.;
    for (size_t ib = 0; ib + 1 < brk.size(); ib++) {
        const REAL PA = brk[ib], PL = brk[ib + 1] - brk[ib];
        if (!(PL > 0.)) continue;
        for (size_t it = 0; it < q.xt.size(); it++) {
            const REAL TH = PA + PL * q.xt[it], WT = PL * q.wt[it];
            const REAL RS = TH < m.thO ? m.Cy / std::sin(TH) : (TH < m.thF ? m.D / std::sin(TH + g.beta) : (H + m.Cy) / std::sin(TH));
            const REAL RP = std::max(m.r0 * std::exp((TH - t1) * k), RS);
            const REAL cT = std::cos(TH), sT = std::sin(TH);
            rb.clear();
            for (int j = 0; j <= q.nRho; j++) rb.push_back(RS + (RP - RS) * (REAL(j) / q.nRho));
            for (const PCircle &c : circles) {
                const REAL dx = m.Cx - c.qx, dy = m.Cy + c.qy;
                const REAL bb = dx * cT + dy * sT, disc = bb * bb - (dx * dx + dy * dy - c.R * c.R);
                const REAL sq = std::sqrt(std::max(disc, 0.));
                for (REAL root : {bb - sq, bb + sq}) rb.push_back(disc > 0. ? std::min(std::max(root, RS), RP) : RS);
            }
            std::sort(rb.begin(), rb.end());
            for (size_t ir = 0; ir + 1 < rb.size(); ir++) {
                const REAL SA = rb[ir], SL = rb[ir + 1] - rb[ir];
                if (!(SL > 0.)) continue;
                for (size_t iq = 0; iq < q.xr.size(); iq++) {
                    const REAL rho = SA + SL * q.xr[iq], W = WT * SL * q.wr[iq];
                    REAL f[2];
                    fpaper(m.Cx - rho * cT, -m.Cy + rho * sT, f);
                    total += (f[0] * sT + f[1] * cT) * rho * rho * W;
                }
            }
        }
    }
    return total;
}

/// P_u / omega = -closed int u U.n dS for f = -grad u (div U = 0): upaper(x, y_paper) (NaN outside the domain of u).
/// Spiral: tan(phi) int u r^2 dtheta; ground surface: crest A-O, face O-F (F = B for I, T for II; panel break at the
/// water line hw if hw > 0), toe ground T-B; nPanels x q Gauss points per piece. Surface points where u is undefined
/// (outside a FE mesh by round-off) are re-evaluated 1e-10 H inside the soil.
template <class U>
REAL BoundaryPower(const Mechanism &m, const Geometry &g, const U &upaper, REAL hw, int nPanels = 24, int q = 8) {
    std::vector<REAL> x01, w01, tq, wq;
    Gauss01(q, x01, w01);
    for (int p = 0; p < nPanels; p++) {
        const REAL a = p == 0 ? 0. : REAL(p) / nPanels, b = p + 1 == nPanels ? 1. : REAL(p + 1) / nPanels;
        for (int i = 0; i < q; i++) tq.push_back(a + (b - a) * x01[i]), wq.push_back((b - a) * w01[i]);
    }
    const REAL eps = 1.e-10 * g.H, H = g.H, sb = g.sb, cb = g.cb;
    auto usurf = [&](REAL x, REAL y, REAL nx, REAL ny) {
        REAL u = upaper(x, y);
        if (!std::isfinite(u)) u = upaper(x - eps * nx, y - eps * ny);
        return u;
    };
    REAL total = 0.;
    if (g.k > 0.) { // spiral (points strictly inside the soil)
        REAL s = 0.;
        for (size_t i = 0; i < tq.size(); i++) {
            const REAL th = m.t1 + (m.t2 - m.t1) * tq[i], r = m.r0 * std::exp((th - m.t1) * g.k);
            s += upaper(m.Cx - r * std::cos(th), -m.Cy + r * std::sin(th)) * r * r * wq[i];
        }
        total += g.k * s * (m.t2 - m.t1);
    }
    { // crest A -> O (y = 0, outward normal (0, -1), U.n = -(Cx - x))
        REAL s = 0.;
        for (size_t i = 0; i < tq.size(); i++) {
            const REAL x = m.Ax * (1. - tq[i]);
            s += usurf(x, 0., 0., -1.) * (-(m.Cx - x)) * wq[i];
        }
        total -= s * m.L;
    }
    { // face O -> F, outward normal (sin b, -cos b)
        const REAL lF = (m.II ? H : m.By) / sb;
        std::vector<REAL> breaks = {0., lF};
        if (hw > 0.) breaks.insert(breaks.begin() + 1, std::min(std::max(hw / sb, 0.), lF));
        for (size_t j = 0; j + 1 < breaks.size(); j++) {
            const REAL a = breaks[j], ll = breaks[j + 1] - a;
            REAL s = 0.;
            for (size_t i = 0; i < tq.size(); i++) {
                const REAL sq = a + ll * tq[i], x = sq * cb, y = sq * sb;
                const REAL un = (y + m.Cy) * sb - (m.Cx - x) * cb;
                s += usurf(x, y, sb, -cb) * un * wq[i];
            }
            total -= s * ll;
        }
    }
    if (m.II) { // toe ground T -> B (y = H)
        const REAL d = m.Bx - g.xtoe;
        REAL s = 0.;
        for (size_t i = 0; i < tq.size(); i++) {
            const REAL x = g.xtoe + d * tq[i];
            s += usurf(x, H, 0., -1.) * (-(m.Cx - x)) * wq[i];
        }
        total -= s * d;
    }
    return total;
}

enum class PuMethod { ENone, EDomain, EBoundary };

inline const char *PuMethodName(PuMethod p) {
    return p == PuMethod::ENone ? "none" : (p == PuMethod::EDomain ? "domain" : "boundary");
}

struct Powers {
    REAL Pmr = 0., Pgamma = 0., Pu = 0.;
};

/// Slope + soil + seepage field: rates of work and Gamma(theta1, theta2, s) (thread safe: const and read-only)
class Problem {
public:
    Geometry geo;
    REAL c = 10., gamma = 20., gammaw = 9.81, gammap = 10.19; ///< gamma' = gamma - gamma_w
    Seepage seep;
    PuMethod pu = PuMethod::ENone;
    REAL dmax = 10.; ///< mechanism II: d <= dmax H

    /// method: "auto" (boundary formula for gradient fields with u, domain quadrature otherwise), "domain", "boundary"
    Problem(REAL betaDeg, REAL H, REAL c_, REAL phiDeg, REAL gamma_, REAL gammaw_, const Seepage &sp,
            const std::string &method = "auto", REAL dmax_ = 10.)
        : geo(H, betaDeg, phiDeg), c(c_), gamma(gamma_), gammaw(gammaw_), gammap(gamma_ - gammaw_), seep(sp), dmax(dmax_) {
        if (seep.None()) pu = PuMethod::ENone;
        else if (method == "domain" || (method == "auto" && !(seep.gradient && seep.u))) pu = PuMethod::EDomain;
        else pu = PuMethod::EBoundary;
        if ((pu == PuMethod::EDomain && !seep.force) || (pu == PuMethod::EBoundary && !seep.u) ||
            (method != "auto" && method != "domain" && method != "boundary")) {
            std::cerr << "la::Problem: P_u method " << method << " not available for this field\n";
            DebugStop();
        }
    }

    /// search box of a mechanism class ('I' or 'II'): theta1, theta2, s
    void Bounds(int kind, REAL lb[3], REAL ub[3]) const {
        const REAL b = geo.betaDeg * M_PI / 180.;
        if (kind == 1) lb[0] = EPS_ANG, lb[1] = EPS_ANG, lb[2] = 1.e-4, ub[0] = M_PI - EPS_ANG, ub[1] = M_PI - b, ub[2] = 1.;
        else lb[0] = EPS_ANG, lb[1] = EPS_ANG, lb[2] = 1., ub[0] = M_PI - EPS_ANG, ub[1] = M_PI - EPS_ANG, ub[2] = 1. + dmax;
    }

    Mechanism Mech(const REAL x[3]) const { return MakeMechanism(geo, x[0], x[1], x[2]); }

    /// f in paper coordinates
    void ForcePaper(REAL x, REAL yp, REAL f[2]) const {
        TPZManVector<REAL, 3> X = {x, -yp, 0.};
        seep.force(X, f);
        f[1] = -f[1];
    }
    /// u in paper coordinates (NaN outside)
    REAL UPaper(REAL x, REAL yp) const {
        TPZManVector<REAL, 3> X = {x, -yp, 0.};
        REAL u;
        return seep.u(X, u) ? u : std::numeric_limits<REAL>::quiet_NaN();
    }

    /// P_u of an admissible mechanism with the method of the problem (or the one given)
    REAL Pu(const Mechanism &m, const Quadrature &q, PuMethod method) const {
        if (method == PuMethod::ENone) return 0.;
        if (method == PuMethod::EBoundary)
            return BoundaryPower(m, geo, [this](REAL x, REAL y) { return UPaper(x, y); }, seep.hw, q.BoundaryPanels(), 8);
        return DomainPower(m, geo, q, [this](REAL x, REAL y, REAL f[2]) { ForcePaper(x, y, f); },
                           FieldCircles(seep, q.nCircleGrade));
    }
    /// P_gamma by the polar quadrature (uniform field (0, gamma') in paper coordinates)
    REAL PgammaQuadrature(const Mechanism &m, const Quadrature &q) const {
        const REAL gp = gammap;
        return DomainPower(m, geo, q, [gp](REAL, REAL, REAL f[2]) { f[0] = 0., f[1] = gp; }, {});
    }
    /// (P_mr, P_gamma, P_u) per unit omega of an admissible mechanism
    Powers Evaluate(const Mechanism &m, const Quadrature &q, bool closedPgamma = true) const {
        Powers p;
        p.Pmr = Pmr(m, geo, c);
        p.Pgamma = closedPgamma ? PgammaClosed(m, geo, gammap) : PgammaQuadrature(m, q);
        p.Pu = Pu(m, q, pu);
        return p;
    }

    /// Gamma = P_mr / (P_gamma + P_u) at x = (theta1, theta2, s); BIG if inadmissible, of the wrong class (kind 1 =
    /// I: s <= 1, 2 = II: s >= 1, 0 = any) or with P_gamma + P_u <= 0 (or not finite)
    REAL GammaFactor(const REAL x[3], int kind, const Quadrature &q) const {
        const Mechanism m = Mech(x);
        bool ok = m.ok;
        if (kind == 1) ok = ok && x[2] <= 1.;
        else if (kind == 2) ok = ok && x[2] >= 1.;
        if (!ok) return BIG;
        const Powers p = Evaluate(m, q);
        const REAL Pext = p.Pgamma + p.Pu;
        const REAL G = Pext > 0. ? p.Pmr / Pext : BIG;
        return std::isfinite(G) ? G : BIG;
    }
    bool Feasible(const REAL x[3], int kind) const {
        return Mech(x).ok && (kind == 1 ? x[2] <= 1. : x[2] >= 1.);
    }
};

/// deterministic uniform random numbers (portable: 53 random bits of mt19937_64)
class Rng {
    std::mt19937_64 fGen;

public:
    explicit Rng(uint64_t seed) : fGen(seed) {}
    REAL Uniform() { return REAL(fGen() >> 11) * 0x1.0p-53; }
    /// uniform integer in [0, n)
    size_t Index(size_t n) {
        const uint64_t lim = UINT64_MAX - UINT64_MAX % n;
        uint64_t v;
        do v = fGen();
        while (v >= lim);
        return size_t(v % n);
    }
};

struct OptPoint {
    REAL x[3] = {0., 0., 0.};
    REAL f = BIG;
};

/// Global-best Particle Swarm Optimisation on a box (limit_analysis.py pso): fun(x) -> value (BIG if inadmissible),
/// feasible(x) -> cheap geometric test used to draw admissible initial particles (up to 50 x 2000 n candidates; the
/// first 4 n admissible ones form the pool, its best n / 2 and n / 2 random others the swarm). Deterministic.
/// Difference from the Python: when fewer than n / 2 mechanisms of the pool have a finite objective (P_ext > 0 is
/// rare for flat slopes), up to extraPools further pools of 4 n admissible draws are evaluated until n / 2 are found;
/// limit_analysis.py keeps the first pool and gives up (Gamma = inf) when none of its mechanisms has P_ext > 0.
inline OptPoint PSO(const std::function<REAL(const REAL *)> &fun, const std::function<bool(const REAL *)> &feasible,
                    const REAL lb[3], const REAL ub[3], int n, int nIter, uint64_t seed, int64_t &nEval,
                    int extraPools = 25, REAL w = 0.7298, REAL c1 = 1.49618, REAL c2 = 1.49618, int stall = 40,
                    REAL vmaxFrac = 0.2) {
    Rng rng(seed);
    const int dim = 3;
    REAL span[3];
    for (int d = 0; d < dim; d++) span[d] = ub[d] - lb[d];
    OptPoint none;
    nEval = 0;
    std::vector<std::array<REAL, 3>> pool;
    std::vector<REAL> fp;
    std::vector<size_t> good;
    const int64_t maxDraw = int64_t(50) * 2000 * n;
    int64_t drawn = 0;
    for (int round = 0; round <= extraPools; round++) {
        const size_t first = pool.size();
        for (; drawn < maxDraw && pool.size() - first < size_t(4 * n); drawn++) {
            std::array<REAL, 3> x;
            for (int d = 0; d < dim; d++) x[d] = lb[d] + rng.Uniform() * span[d];
            if (feasible(x.data())) pool.push_back(x);
        }
        for (size_t i = first; i < pool.size(); i++) {
            fp.push_back(fun(pool[i].data()));
            if (fp[i] < BIG) good.push_back(i);
        }
        nEval += pool.size() - first;
        if (good.size() >= size_t(n / 2) || pool.size() == first || drawn >= maxDraw) break;
    }
    if (pool.empty()) return none;
    if (good.empty()) {
        std::copy(pool[0].begin(), pool[0].end(), none.x);
        return none;
    }
    std::stable_sort(good.begin(), good.end(), [&fp](size_t a, size_t b) { return fp[a] < fp[b]; });
    const size_t nb = std::min(size_t(n / 2), good.size());
    std::vector<size_t> pick(good.begin(), good.begin() + nb), rest(good.begin() + nb, good.end());
    const size_t nr = std::min(rest.size(), size_t(n) - nb);
    for (size_t i = 0; i < nr; i++) { // random subset without replacement (partial Fisher-Yates)
        const size_t j = i + rng.Index(rest.size() - i);
        std::swap(rest[i], rest[j]);
        pick.push_back(rest[i]);
    }
    std::vector<std::array<REAL, 3>> X, V(n), P;
    std::vector<REAL> F, PF;
    for (size_t i : pick) X.push_back(pool[i]), F.push_back(fp[i]);
    while (int(X.size()) < n) { // fill with random points of the box
        std::array<REAL, 3> x;
        for (int d = 0; d < dim; d++) x[d] = lb[d] + rng.Uniform() * span[d];
        X.push_back(x), F.push_back(fun(x.data())), nEval++;
    }
    for (int i = 0; i < n; i++)
        for (int d = 0; d < dim; d++) V[i][d] = (rng.Uniform() - 0.5) * 0.2 * span[d];
    P = X, PF = F;
    int g = int(std::min_element(PF.begin(), PF.end()) - PF.begin());
    REAL bestLast = PF[g];
    int itLast = 0;
    for (int it = 0; it < nIter; it++) {
        for (int i = 0; i < n; i++)
            for (int d = 0; d < dim; d++) {
                const REAL r1 = rng.Uniform(), r2 = rng.Uniform(), vmax = vmaxFrac * span[d];
                V[i][d] = w * V[i][d] + c1 * r1 * (P[i][d] - X[i][d]) + c2 * r2 * (P[g][d] - X[i][d]);
                V[i][d] = std::min(std::max(V[i][d], -vmax), vmax);
                X[i][d] = std::min(std::max(X[i][d] + V[i][d], lb[d]), ub[d]);
            }
        for (int i = 0; i < n; i++) {
            F[i] = fun(X[i].data());
            if (F[i] < PF[i]) P[i] = X[i], PF[i] = F[i];
        }
        nEval += n;
        g = int(std::min_element(PF.begin(), PF.end()) - PF.begin());
        if (PF[g] < bestLast * (1. - 1.e-9)) itLast = it;
        bestLast = PF[g];
        if (it - itLast > stall) break;
    }
    OptPoint r;
    std::copy(P[g].begin(), P[g].end(), r.x);
    r.f = PF[g];
    return r;
}

/// Bounded Nelder-Mead of scipy.optimize (minimize(method="Nelder-Mead", bounds=..., initial_simplex=...)): standard
/// coefficients, trial points clipped to the box, stop when max|x_i - x_0| <= xatol and max|f_i - f_0| <= fatol
inline OptPoint NelderMeadScipy(const std::function<REAL(const REAL *)> &fun, const REAL x0[3], REAL sim0[4][3],
                                const REAL lb[3], const REAL ub[3], REAL xatol, REAL fatol, int maxiter, int maxfev,
                                int64_t &nEval) {
    (void)x0; // the simplex sim0 carries the start point (row 0)
    const int N = 3;
    const REAL rho = 1., chi = 2., psi = 0.5, sigma = 0.5;
    std::array<std::array<REAL, 3>, 4> sim;
    std::array<REAL, 4> fsim;
    auto clip = [&](std::array<REAL, 3> &x) {
        for (int d = 0; d < N; d++) x[d] = std::min(std::max(x[d], lb[d]), ub[d]);
    };
    for (int j = 0; j <= N; j++)
        for (int d = 0; d < N; d++) {
            REAL v = sim0[j][d];
            if (v > ub[d]) v = 2. * ub[d] - v; // reflect into the interior, then clip
            sim[j][d] = std::min(std::max(v, lb[d]), ub[d]);
        }
    int fcalls = 0;
    auto f = [&](const std::array<REAL, 3> &x) {
        fcalls++;
        return fun(x.data());
    };
    auto sortSim = [&]() {
        std::array<int, 4> idx = {0, 1, 2, 3};
        std::stable_sort(idx.begin(), idx.end(), [&fsim](int a, int b) { return fsim[a] < fsim[b]; });
        auto s2 = sim;
        auto f2 = fsim;
        for (int j = 0; j <= N; j++) sim[j] = s2[idx[j]], fsim[j] = f2[idx[j]];
    };
    for (int j = 0; j <= N; j++) fsim[j] = f(sim[j]);
    sortSim();
    int iterations = 1;
    while (fcalls < maxfev && iterations < maxiter) {
        REAL dx = 0., df = 0.;
        for (int j = 1; j <= N; j++) {
            for (int d = 0; d < N; d++) dx = std::max(dx, std::fabs(sim[j][d] - sim[0][d]));
            df = std::max(df, std::fabs(fsim[0] - fsim[j]));
        }
        if (dx <= xatol && df <= fatol) break;
        std::array<REAL, 3> xbar = {0., 0., 0.}, xr, xe, xc, xcc;
        for (int j = 0; j < N; j++)
            for (int d = 0; d < N; d++) xbar[d] += sim[j][d];
        for (int d = 0; d < N; d++) xbar[d] /= N;
        for (int d = 0; d < N; d++) xr[d] = (1. + rho) * xbar[d] - rho * sim[N][d];
        clip(xr);
        const REAL fxr = f(xr);
        bool shrink = false;
        if (fxr < fsim[0]) {
            for (int d = 0; d < N; d++) xe[d] = (1. + rho * chi) * xbar[d] - rho * chi * sim[N][d];
            clip(xe);
            const REAL fxe = f(xe);
            if (fxe < fxr) sim[N] = xe, fsim[N] = fxe;
            else sim[N] = xr, fsim[N] = fxr;
        } else if (fxr < fsim[N - 1]) {
            sim[N] = xr, fsim[N] = fxr;
        } else if (fxr < fsim[N]) { // outside contraction
            for (int d = 0; d < N; d++) xc[d] = (1. + psi * rho) * xbar[d] - psi * rho * sim[N][d];
            clip(xc);
            const REAL fxc = f(xc);
            if (fxc <= fxr) sim[N] = xc, fsim[N] = fxc;
            else shrink = true;
        } else { // inside contraction
            for (int d = 0; d < N; d++) xcc[d] = (1. - psi) * xbar[d] + psi * sim[N][d];
            clip(xcc);
            const REAL fxcc = f(xcc);
            if (fxcc < fsim[N]) sim[N] = xcc, fsim[N] = fxcc;
            else shrink = true;
        }
        if (shrink)
            for (int j = 1; j <= N; j++) {
                for (int d = 0; d < N; d++) sim[j][d] = sim[0][d] + sigma * (sim[j][d] - sim[0][d]);
                clip(sim[j]);
                fsim[j] = f(sim[j]);
            }
        iterations++;
        sortSim();
    }
    nEval += fcalls;
    OptPoint r;
    std::copy(sim[0].begin(), sim[0].end(), r.x);
    r.f = *std::min_element(fsim.begin(), fsim.end());
    return r;
}

/// restarted bounded Nelder-Mead polish of limit_analysis.py (nelder_mead): simplex steps (0.02, 0.02, 0.01), then
/// x 0.2 per restart; a restart is accepted if it does not increase the value
inline OptPoint Polish(const std::function<REAL(const REAL *)> &fun, const REAL x0[3], const REAL lb[3], const REAL ub[3],
                       int64_t &nEval, int restarts = 2) {
    OptPoint best;
    std::copy(x0, x0 + 3, best.x);
    for (int d = 0; d < 3; d++) best.x[d] = std::min(std::max(best.x[d], lb[d]), ub[d]);
    best.f = fun(best.x);
    nEval++;
    REAL step[3] = {0.02, 0.02, 0.01};
    for (int r = 0; r < restarts; r++) {
        REAL sim[4][3];
        for (int d = 0; d < 3; d++) sim[0][d] = best.x[d];
        for (int i = 0; i < 3; i++) {
            for (int d = 0; d < 3; d++) sim[i + 1][d] = best.x[d];
            sim[i + 1][i] = best.x[i] + step[i] <= ub[i] ? best.x[i] + step[i] : best.x[i] - step[i];
        }
        const OptPoint res = NelderMeadScipy(fun, best.x, sim, lb, ub, 1.e-10, 1.e-13, 4000, 8000, nEval);
        if (res.f <= best.f) best = res;
        for (REAL &s : step) s *= 0.2;
    }
    return best;
}

struct Settings {
    std::vector<int> kinds = {1, 2};                 ///< mechanism classes searched: 1 = I, 2 = II
    std::vector<uint64_t> seeds = {0, 1, 2};         ///< one PSO run per class and seed
    int nParticles = 40, nIter = 150;
    std::string quadSearch = "coarse", quadFinal = "fine";
    bool polish = true;
    int extraPools = 25;                             ///< further initial pools when P_ext > 0 is rare (0: as the Python)
    int nThreads = 0;                                ///< 0: hardware concurrency
};

struct Run {
    int kind = 1;
    uint64_t seed = 0;
    bool found = false;
    REAL xPSO[3] = {0., 0., 0.}, x[3] = {0., 0., 0.};
    REAL GammaPSO = BIG, Gamma = BIG;
    int64_t nEval = 0;
    double seconds = 0.;
};

struct Result {
    bool found = false;                       ///< false: Gamma = inf (no admissible mechanism with P_ext > 0)
    REAL Gamma = INFINITY, Hcrit = INFINITY;
    int kind = 0;                             ///< 1 = I (B on the face), 2 = II (B on the toe ground)
    REAL x[3] = {0., 0., 0.};                 ///< theta1, theta2, s (eta = s for I, d / H = s - 1 for II)
    Mechanism mech;
    Powers P;                                 ///< per unit omega, quadrature quadFinal
    REAL A[2] = {0., 0.}, B[2] = {0., 0.}, C[2] = {0., 0.}; ///< paper coordinates (y down); NeoPZ: (x, -y)
    std::vector<Run> runs;
    std::vector<std::string> atBound;         ///< parameters at a search bound
    bool dmaxReached = false;
    REAL seedSpread = 0.;                     ///< max / min - 1 of the runs of the best class
    REAL bestByClass[3] = {INFINITY, INFINITY, INFINITY};
    int64_t nEval = 0;
    double seconds = 0.;
    PuMethod pu = PuMethod::ENone;
};

/// Upper-bound stability factor Gamma = min P_mr / (P_gamma + P_u) (Eq. 55); H_crit = Gamma H
inline Result StabilityFactor(const Problem &prob, const Settings &st = Settings()) {
    const auto t0 = std::chrono::steady_clock::now();
    const Quadrature qs = Quadrature::Level(st.quadSearch), qf = Quadrature::Level(st.quadFinal);
    if (!qs.Valid() || !qf.Valid()) {
        std::cerr << "la::StabilityFactor: unknown quadrature level " << st.quadSearch << " / " << st.quadFinal << "\n";
        DebugStop();
    }
    Result res;
    res.pu = prob.pu;
    for (int kind : st.kinds)
        for (uint64_t seed : st.seeds) {
            Run r;
            r.kind = kind, r.seed = seed;
            res.runs.push_back(r);
        }
    auto work = [&](Run &r) {
        const auto ts = std::chrono::steady_clock::now();
        REAL lb[3], ub[3];
        prob.Bounds(r.kind, lb, ub);
        const int kind = r.kind;
        const auto fs = [&prob, &qs, kind](const REAL *x) { return prob.GammaFactor(x, kind, qs); };
        const auto ff = [&prob, &qf, kind](const REAL *x) { return prob.GammaFactor(x, kind, qf); };
        const auto feas = [&prob, kind](const REAL *x) { return prob.Feasible(x, kind); };
        const OptPoint p = PSO(fs, feas, lb, ub, st.nParticles, st.nIter, r.seed, r.nEval, st.extraPools);
        if (p.f < BIG) {
            r.found = true;
            std::copy(p.x, p.x + 3, r.xPSO);
            r.GammaPSO = p.f;
            if (st.polish) {
                const OptPoint q = Polish(ff, p.x, lb, ub, r.nEval);
                std::copy(q.x, q.x + 3, r.x);
                r.Gamma = q.f;
            } else {
                std::copy(p.x, p.x + 3, r.x);
                r.Gamma = ff(p.x);
                r.nEval++;
            }
            if (!(r.Gamma < BIG)) r.found = false;
        }
        r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts).count();
    };
    int nth = st.nThreads > 0 ? st.nThreads : int(std::max(1u, std::thread::hardware_concurrency()));
    nth = std::min<int>(nth, int(res.runs.size()));
    if (nth <= 1) {
        for (Run &r : res.runs) work(r);
    } else {
        std::atomic<size_t> next(0);
        std::vector<std::thread> th;
        for (int t = 0; t < nth; t++)
            th.emplace_back([&]() {
                for (size_t i = next++; i < res.runs.size(); i = next++) work(res.runs[i]);
            });
        for (auto &t : th) t.join();
    }
    const Run *best = nullptr;
    for (const Run &r : res.runs) {
        res.nEval += r.nEval;
        if (r.found && (!best || r.Gamma < best->Gamma)) best = &r;
        if (r.found) res.bestByClass[r.kind] = std::min(res.bestByClass[r.kind], r.Gamma);
    }
    res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!best) return res;
    res.found = true;
    std::copy(best->x, best->x + 3, res.x);
    res.mech = prob.Mech(res.x);
    res.P = prob.Evaluate(res.mech, qf);
    res.Gamma = res.P.Pmr / (res.P.Pgamma + res.P.Pu);
    res.Hcrit = res.Gamma * prob.geo.H;
    res.kind = res.x[2] > 1. ? 2 : 1;
    const Mechanism &m = res.mech;
    res.A[0] = m.Ax, res.A[1] = 0.;
    res.B[0] = m.Bx, res.B[1] = m.By;
    res.C[0] = m.Cx, res.C[1] = -m.Cy;
    REAL lb[3], ub[3];
    prob.Bounds(best->kind, lb, ub);
    static const char *names[3] = {"theta1", "theta2", "s"};
    for (int i = 0; i < 3; i++) {
        const REAL tol = 1.e-6 * (ub[i] - lb[i]);
        if (std::min(std::fabs(res.x[i] - lb[i]), std::fabs(res.x[i] - ub[i])) < tol)
            res.atBound.push_back(std::string(names[i]) + (std::fabs(res.x[i] - lb[i]) < tol ? "=lower" : "=upper"));
    }
    res.dmaxReached = res.kind == 2 && res.x[2] > 1. + 0.999 * prob.dmax;
    REAL gmin = INFINITY, gmax = 0.;
    for (const Run &r : res.runs)
        if (r.found && r.kind == best->kind) gmin = std::min(gmin, r.Gamma), gmax = std::max(gmax, r.Gamma);
    res.seedSpread = gmax / gmin - 1.;
    return res;
}

} // namespace la
} // namespace slope

#endif
