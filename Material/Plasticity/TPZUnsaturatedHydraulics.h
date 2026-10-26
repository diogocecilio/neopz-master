/**
 * @file TPZUnsaturatedHydraulics.h
 * @brief Constitutive relations of an unsaturated soil with water and air: retention curve S_w(p_c) and relative
 * permeabilities k_rw(S_w), k_ra(S_w), used by TPZMatUnsaturatedUPwPa.
 */

#ifndef TPZUNSATURATEDHYDRAULICS_H
#define TPZUNSATURATEDHYDRAULICS_H

#include "pzreal.h"
#include <algorithm>
#include <cmath>
#include <iostream>

/**
 * @ingroup material
 * @brief Retention curve and relative permeabilities of the two fluid phases (Sects. 3 and 5 of the notes
 * "Formulação poroelástica não saturada")
 *
 * The capillary pressure is \f$p_c = p_a - p_w\f$ (kPa); for \f$p_c \le 0\f$ the soil is saturated, \f$S_w = 1\f$.
 *
 * Retention curves:
 *  - EFredlundXing: \f$S_w = \left[\ln\left(e + (p_c/a)^n\right)\right]^{-m}\f$ with the correction factor
 *    \f$C(\psi) = 1\f$ (Fredlund and Xing 1994, as recommended by Leong and Rahardjo 1997 and used by
 *    Rahardjo et al. 2007); \f$S_w = \theta_w/\theta_s\f$;
 *  - EVanGenuchten: \f$S_e = \left[1 + (a_{VG} p_c)^{b_{VG}}\right]^{-m_{VG}}\f$, \f$m_{VG} = 1 - 1/b_{VG}\f$,
 *    \f$S_w = S_{wr} + (1 - S_{wr}) S_e\f$.
 *
 * Relative permeability of the water:
 *  - EPowerLaw: \f$k_{rw} = S_e^{p}\f$ (Leong and Rahardjo 1997b, \f$k_w = k_s\Theta^p\f$);
 *  - EMualem: \f$k_{rw} = S_e^{1/2}\left\{1 - \left[1 - S_e^{1/m_{VG}}\right]^{m_{VG}}\right\}^2\f$.
 *
 * Relative permeability of the air (non-wetting phase):
 *  - EBrooksCorey: \f$k_{ra} = (1 - S_e)^2\left(1 - S_e^{(2+\lambda)/\lambda}\right)\f$;
 *  - EMualemAir: \f$k_{ra} = (1 - S_e)^{1/2}\left[1 - S_e^{1/m_{VG}}\right]^{2m_{VG}}\f$.
 *
 * The relative permeabilities are bounded from below by fKrwMin and fKraMin, which keep the flow blocks
 * of the tangent regular where a phase disappears (Sect. 13.9 of the notes).
 */
struct TPZUnsaturatedHydraulics {
    enum ERetention { EFredlundXing = 0, EVanGenuchten = 1 };
    enum EWaterPermeability { EPowerLaw = 0, EMualem = 1 };
    enum EAirPermeability { EBrooksCorey = 0, EMualemAir = 1 };

    ERetention fRetention = EFredlundXing;
    EWaterPermeability fKrw = EPowerLaw;
    EAirPermeability fKra = EBrooksCorey;

    /** @brief Fredlund-Xing: a (kPa), n, m */
    REAL fA = 50., fN = 1., fM = 1.;
    /** @brief van Genuchten: a_VG (1/kPa), b_VG and residual saturation S_wr */
    REAL fAVG = 0.1, fBVG = 2., fSwr = 0.;
    /** @brief exponent p of the power law of the water */
    REAL fPower = 4.;
    /** @brief exponent lambda of Brooks-Corey (air) */
    REAL fLambdaBC = 3.;
    /** @brief lower bounds of k_rw and k_ra */
    REAL fKrwMin = 1.e-12, fKraMin = 1.e-4;

    /** @brief State of the hydraulic relations at a capillary pressure */
    struct TState {
        REAL S = 1.;    ///< water saturation S_w
        REAL Gamma = 0.; ///< \f$\Gamma = \partial S_w/\partial p_w = -\partial S_w/\partial p_a = -dS_w/dp_c \ge 0\f$
        REAL krw = 1.;  ///< relative permeability of the water
        REAL dkrw = 0.; ///< \f$dk_{rw}/dS_w\f$
        REAL kra = 0.;  ///< relative permeability of the air
        REAL dkra = 0.; ///< \f$dk_{ra}/dS_w\f$
    };

    REAL MVG() const { return 1. - 1. / fBVG; }

    /** @brief Saturation and its derivative \f$\Gamma = -dS_w/dp_c\f$ */
    void Retention(REAL pc, REAL &S, REAL &Gamma) const {
        if (pc <= 0.) {
            S = 1.;
            Gamma = 0.;
            return;
        }
        if (fRetention == EFredlundXing) {
            const REAL r = pc / fA, x = std::pow(r, fN), L = std::log(M_E + x);
            S = std::pow(L, -fM);
            Gamma = fM * std::pow(L, -fM - 1.) * fN * std::pow(r, fN - 1.) / (fA * (M_E + x));
        } else {
            const REAL m = MVG(), y = std::pow(fAVG * pc, fBVG);
            const REAL Se = std::pow(1. + y, -m);
            S = fSwr + (1. - fSwr) * Se;
            // dSe/dpc = -m b a^b pc^(b-1) (1 + y)^(-m-1)
            Gamma = (1. - fSwr) * m * fBVG * std::pow(fAVG, fBVG) * std::pow(pc, fBVG - 1.) * std::pow(1. + y, -m - 1.);
        }
    }

    /** @brief Effective saturation \f$S_e\f$ and \f$dS_e/dS_w\f$ */
    void Effective(REAL S, REAL &Se, REAL &dSe) const {
        const REAL swr = (fRetention == EVanGenuchten) ? fSwr : 0.;
        dSe = 1. / (1. - swr);
        Se = std::min(std::max((S - swr) * dSe, REAL(0.)), REAL(1.));
    }

    /** @brief Relative permeability of the water and its derivative with respect to S_w */
    void WaterPermeability(REAL S, REAL &kr, REAL &dkr) const {
        REAL Se, dSe;
        Effective(S, Se, dSe);
        REAL raw, draw;
        if (fKrw == EPowerLaw) {
            raw = std::pow(Se, fPower);
            draw = fPower * std::pow(Se, fPower - 1.);
        } else {
            const REAL m = MVG();
            const REAL Sec = std::min(Se, REAL(1. - 1.e-10));
            const REAL q = std::pow(Sec, 1. / m), f = 1. - std::pow(1. - q, m);
            raw = std::sqrt(Se) * f * f;
            // df/dSe = m (1-q)^(m-1) q / (m Se) = (1-q)^(m-1) q / Se
            const REAL df = std::pow(1. - q, m - 1.) * q / std::max(Sec, REAL(1.e-300));
            draw = (Se > 0. ? 0.5 / std::sqrt(Se) * f * f : 0.) + std::sqrt(Se) * 2. * f * df;
        }
        if (raw > fKrwMin) {
            kr = raw;
            dkr = draw * dSe;
        } else {
            kr = fKrwMin;
            dkr = 0.;
        }
    }

    /** @brief Relative permeability of the air and its derivative with respect to S_w */
    void AirPermeability(REAL S, REAL &kr, REAL &dkr) const {
        REAL Se, dSe;
        Effective(S, Se, dSe);
        REAL raw, draw;
        if (fKra == EBrooksCorey) {
            const REAL e = (2. + fLambdaBC) / fLambdaBC;
            raw = (1. - Se) * (1. - Se) * (1. - std::pow(Se, e));
            draw = -2. * (1. - Se) * (1. - std::pow(Se, e)) - (1. - Se) * (1. - Se) * e * std::pow(Se, e - 1.);
        } else {
            const REAL m = MVG();
            const REAL Sec = std::max(Se, REAL(1.e-12));
            const REAL q = std::pow(Sec, 1. / m), g = std::pow(1. - q, 2. * m);
            raw = std::sqrt(1. - Se) * g;
            const REAL dg = -2. * std::pow(1. - q, 2. * m - 1.) * q / Sec;
            draw = (Se < 1. ? -0.5 / std::sqrt(1. - Se) * g : 0.) + std::sqrt(1. - Se) * dg;
        }
        if (raw > fKraMin) {
            kr = raw;
            dkr = draw * dSe;
        } else {
            kr = fKraMin;
            dkr = 0.;
        }
    }

    /** @brief All the relations at the pressures p_w, p_a */
    TState Evaluate(REAL pw, REAL pa) const {
        TState st;
        Retention(pa - pw, st.S, st.Gamma);
        WaterPermeability(st.S, st.krw, st.dkrw);
        AirPermeability(st.S, st.kra, st.dkra);
        return st;
    }

    void Print(std::ostream &out) const {
        out << "  retention " << (fRetention == EFredlundXing ? "Fredlund-Xing" : "van Genuchten");
        if (fRetention == EFredlundXing) out << " a = " << fA << " n = " << fN << " m = " << fM;
        else out << " a_VG = " << fAVG << " b_VG = " << fBVG << " S_wr = " << fSwr;
        out << "\n  k_rw " << (fKrw == EPowerLaw ? "power law p = " : "Mualem") << (fKrw == EPowerLaw ? fPower : 0.)
            << " (min " << fKrwMin << "), k_ra " << (fKra == EBrooksCorey ? "Brooks-Corey lambda = " : "Mualem")
            << (fKra == EBrooksCorey ? fLambdaBC : 0.) << " (min " << fKraMin << ")\n";
    }
};

#endif
