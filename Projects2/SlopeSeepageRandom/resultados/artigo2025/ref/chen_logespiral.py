#!/usr/bin/env python3
"""Limite superior log-espiral de Chen (1975), mecanismo pelo pé, talude simples com crista horizontal.

Ns = γ H / c mínimo sobre (θ0, θh); sem percolação. Serve de verificação independente dos valores secos do artigo
(Vargas Ceron et al. 2025): Γ = 1,777 do exemplo c-φ de Cho e H_crit da Fig. 8 em h_w/H = 0.
"""
import numpy as np
from scipy.optimize import minimize


def ns(beta_deg, phi_deg):
    b = np.radians(beta_deg)
    t = np.tan(np.radians(phi_deg))

    def f(x):
        t0, th = x
        if not 0 < t0 < th < np.pi:
            return 1e9
        e = np.exp((th - t0) * t)
        hr = np.sin(th) * e - np.sin(t0)  # H / r0
        lr = np.sin(th - t0) / np.sin(th) - np.sin(th + b) / (np.sin(th) * np.sin(b)) * hr  # L / r0
        if hr <= 0 or lr <= 0:
            return 1e9
        f1 = ((3 * t * np.cos(th) + np.sin(th)) * e ** 3 - 3 * t * np.cos(t0) - np.sin(t0)) / (3 * (1 + 9 * t * t))
        f2 = lr * (2 * np.cos(t0) - lr) * np.sin(t0) / 6
        f3 = e * (np.sin(th - t0) - lr * np.sin(th)) * (np.cos(t0) - lr + np.cos(th) * e) / 6
        w = f1 - f2 - f3
        if w <= 0:
            return 1e9
        return (e * e - 1) / (2 * t) / w * hr

    best, x0 = 1e9, None
    for t0 in np.linspace(0.02, 2.5, 60):
        for th in np.linspace(t0 + 0.02, 3.1, 60):
            v = f((t0, th))
            if v < best:
                best, x0 = v, (t0, th)
    if x0 is None:
        return float("inf")
    r = minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-10, "maxiter": 20000})
    return r.fun if r.fun < 1e8 else float("inf")


if __name__ == "__main__":
    gw = 9.81
    print(f"Cho c-φ (β = 45°, φ = 30°, c = 10, γ = 20, H = 10 m): Γ = {ns(45, 30) * 10 / (20 * 10):.4f}")
    for nome, c, phi in (("A (c = 6, φ = 32°)", 6., 32.), ("B (c = 11,7, φ = 24,7°)", 11.7, 24.7)):
        for b in (30, 35, 60):
            n = ns(b, phi)
            print(f"solo {nome}, β = {b}°: Ns = {n:.3f}, H_crit = {n * c / (18 - gw):.2f} m (γ' = {18 - gw:.2f})")
