"""Os eixos das figuras de Vargas Ceron et al. (2025), "Stability Analysis of Slope Subjected to Seepage Forces Considering
Spatial Variability of Soil Properties", Int J Numer Anal Methods Geomech 49:2459-2491 (doi 10.1002/nag.3993), lidos nas
figuras do artigo: limites, passo das marcacoes, escala (linear ou log), rotulos e estilos das curvas. As figuras do
relatorio usam estes eixos; quando os dados nao cabem, amplia() estende o eixo so o necessario para conter a faixa do
Ceron e a dos dados, e a legenda da figura diz qual era a faixa original.

Cada painel: {"x": (min, max, passo), "y": (min, max, passo), "xlog": bool, "xrot": rotulo, "yrot": rotulo}.
Rotulos em portugues e em LaTeX (mathtext); onde o Ceron escreve Gamma (fator de estabilidade da analise limite), o
relatorio escreve FS (o fator de seguranca do SRM)."""

import math

# estilos das curvas das PDFs, na ordem em que o Ceron os usa (Figs. 27, 31a, 32a e 33a): cor, tracejado
ESTILOS = [("k", "-"), ("k", "--"), ("r", "-"), ("r", "--"), ("b", "-"), ("b", "--"), ("0.35", "--")]

# os valores de cada familia nas figuras do Ceron, na ordem das legendas: o estilo de um valor e o que o Ceron usa para
# ele (o mesmo valor, a mesma linha); um valor que o Ceron nao tem recebe o proximo estilo livre
VALORES_CERON = {
    "covk": [0.0, 0.6, 0.75, 0.9, 1.0],
    "covc": [0.1, 0.3, 0.5, 0.7],
    "covphi": [0.05, 0.1, 0.15, 0.2],
    "s": [1.0, 1.5, 2.0, 5.0, 10.0, 20.0],
    "alfa": [1.0, 2.0, 3.0, 4.0, 5.0],
    "hw": [0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
}

# cores das barras da Fig. 26 (k_v, phi, c), de cima para baixo
CORES_SPEARMAN = {"k": "#D1344B", "phi": "#2E7591", "c": "#E9A93B"}

PDF = "Densidade de probabilidade"
FIG = {
    # Fig. 25: configuracao de referencia
    "25a": {"x": (0.0, 3.5, 0.5), "y": (0.0, 1.6, 0.2), "xlog": False, "xrot": r"$FS$", "yrot": PDF},
    "25b": {"x": (1.0, 1.0e4, None), "y": (0.0, 1.0, 0.2), "xlog": True, "xrot": "Número de simulações",
            "yrot": r"$P_f$"},
    # Fig. 26: Spearman parcial de 2a ordem (barras horizontais: k_v, phi, c de cima para baixo)
    "26": {"x": (-1.0, 1.0, 0.25), "y": None, "xlog": False, "xrot": r"$\rho_{X\,FS,YZ}$",
           "yrot": "Variáveis aleatórias de entrada"},
    # Fig. 27: PDFs com o CoV de k_v (a), c (b) e phi (c)
    "27a": {"x": (0.5, 3.0, 0.5), "y": (0.0, 1.5, 0.25), "xlog": False, "xrot": r"$FS$", "yrot": PDF,
            "legenda": r"$CoV(k_v)$"},
    "27b": {"x": (0.0, 3.0, 0.5), "y": (0.0, 2.5, 0.5), "xlog": False, "xrot": r"$FS$", "yrot": PDF,
            "legenda": r"$CoV(c)$"},
    "27c": {"x": (0.0, 3.0, 0.5), "y": (0.0, 1.5, 0.25), "xlog": False, "xrot": r"$FS$", "yrot": PDF,
            "legenda": r"$CoV(\varphi)$"},
    # Fig. 28: P_f contra o CoV; em (a), CoV(k_v) = 0 e uma linha tracejada horizontal
    "28a": {"x": (0.6, 1.0, 0.1), "y": (0.110, 0.120, 0.002), "xlog": False, "xrot": r"$CoV(k_v)$", "yrot": r"$P_f$"},
    "28b": {"x": (0.1, 0.7, 0.1), "y": (0.0, 0.4, 0.1), "xlog": False, "xrot": r"$CoV(c)$", "yrot": r"$P_f$"},
    "28c": {"x": (0.05, 0.20, 0.025), "y": (0.10, 0.20, 0.02), "xlog": False, "xrot": r"$CoV(\varphi)$",
            "yrot": r"$P_f$"},
    # Fig. 29: media do fator contra o CoV
    "29a": {"x": (0.6, 1.0, 0.1), "y": (1.1, 1.5, 0.1), "xlog": False, "xrot": r"$CoV(k_v)$", "yrot": r"$\mu(FS)$"},
    "29b": {"x": (0.1, 0.7, 0.1), "y": (1.1, 1.5, 0.1), "xlog": False, "xrot": r"$CoV(c)$", "yrot": r"$\mu(FS)$"},
    "29c": {"x": (0.05, 0.20, 0.025), "y": (1.1, 1.5, 0.1), "xlog": False, "xrot": r"$CoV(\varphi)$",
            "yrot": r"$\mu(FS)$"},
    # Fig. 30: CoV do fator contra o CoV
    "30a": {"x": (0.6, 1.0, 0.1), "y": (0.1, 0.6, 0.1), "xlog": False, "xrot": r"$CoV(k_v)$", "yrot": r"$CoV(FS)$"},
    "30b": {"x": (0.1, 0.7, 0.1), "y": (0.1, 0.6, 0.1), "xlog": False, "xrot": r"$CoV(c)$", "yrot": r"$CoV(FS)$"},
    "30c": {"x": (0.05, 0.20, 0.025), "y": (0.1, 0.6, 0.1), "xlog": False, "xrot": r"$CoV(\varphi)$",
            "yrot": r"$CoV(FS)$"},
    # Fig. 31: distancias de autocorrelacao (L_x, L_y) = s (20 m, 2 m); s = 400 e uma linha tracejada em (b) a (d)
    "31a": {"x": (0.0, 3.0, 0.5), "y": (0.0, 1.5, 0.25), "xlog": False, "xrot": r"$FS$", "yrot": PDF,
            "legenda": r"$s$"},
    "31b": {"x": (0.0, 20.0, 5.0), "y": (0.0, 0.4, 0.1), "xlog": False, "xrot": r"$s$", "yrot": r"$P_f$"},
    "31c": {"x": (0.0, 20.0, 5.0), "y": (0.1, 0.6, 0.1), "xlog": False, "xrot": r"$s$", "yrot": r"$CoV(FS)$"},
    "31d": {"x": (0.0, 20.0, 5.0), "y": (1.1, 1.6, 0.1), "xlog": False, "xrot": r"$s$", "yrot": r"$\mu(FS)$"},
    # Fig. 32: anisotropia alfa = k_h/k_v
    "32a": {"x": (0.0, 4.0, 0.5), "y": (0.0, 1.5, 0.25), "xlog": False, "xrot": r"$FS$", "yrot": PDF,
            "legenda": r"$\alpha$"},
    "32b": {"x": (1.0, 5.0, 1.0), "y": (0.0, 0.150, 0.025), "xlog": False, "xrot": r"$\alpha$", "yrot": r"$P_f$"},
    "32c": {"x": (1.0, 5.0, 1.0), "y": (0.20, 0.30, 0.02), "xlog": False, "xrot": r"$\alpha$",
            "yrot": r"$CoV(FS)$"},
    "32d": {"x": (1.0, 5.0, 1.0), "y": (1.4, 2.0, 0.2), "xlog": False, "xrot": r"$\alpha$", "yrot": r"$\mu(FS)$"},
    # Fig. 33: rebaixamento h_w/H (no relatorio, a intensidade da chuva ocupa o lugar de h_w/H; o eixo x e outro)
    "33a": {"x": (0.5, 3.5, 0.5), "y": (0.0, 1.5, 0.25), "xlog": False, "xrot": r"$FS$", "yrot": PDF,
            "legenda": r"$h_w/H$"},
    "33b": {"x": (0.5, 1.0, 0.1), "y": (0.0, 0.150, 0.025), "xlog": False, "xrot": r"$h_w/H$", "yrot": r"$P_f$"},
    "33c": {"x": (0.5, 1.0, 0.1), "y": (0.20, 0.30, 0.02), "xlog": False, "xrot": r"$h_w/H$",
            "yrot": r"$CoV(FS)$"},
    "33d": {"x": (0.5, 1.0, 0.1), "y": (1.3, 1.7, 0.1), "xlog": False, "xrot": r"$h_w/H$", "yrot": r"$\mu(FS)$"},
}

# Figs. 12 a 14: barras de cor dos campos (horizontais, abaixo de cada painel), escala de cores azul -> ciano -> verde
# -> amarelo -> vermelho -> magenta. A permeabilidade do Ceron esta em m2/(kPa s) (= m4/(kN s), Tabela 2): k (m/s) / gamma_w
CAMPOS = {
    "k": {"lim": (1.0e-8, 4.0e-6), "marcas": [1.0e-8, 1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6],
          "rot": r"Permeabilidade vertical (m$^2$ kPa$^{-1}$ s$^{-1}$)"},
    "c": {"lim": (5.0, 15.0), "marcas": [5.0, 7.5, 10.0, 12.5, 15.0], "rot": r"Coesão (kPa)"},
    "phi": {"lim": (25.0, 40.0), "marcas": [25.0, 27.5, 30.0, 32.5, 35.0, 37.5, 40.0], "rot": r"Ângulo de atrito ($^\circ$)"},
    "u": {"lim": (-50.0, 0.0), "marcas": [-50.0, -40.0, -30.0, -20.0, -10.0, 0.0],
          "rot": r"Excesso de poropressão (kPa)"},
}
CORES_CAMPOS = ["#0000b0", "#0060ff", "#00d0ff", "#20e080", "#a0ff20", "#ffff00", "#ffa000", "#ff2000", "#d00060",
                "#a00090"]

GAMMA_W = 9.81  # kN/m3, para converter k (m/s) na unidade do Ceron


def _passo_bom(passo, faixa, max_marcas=11):
    """o menor multiplo 'redondo' do passo do Ceron (x1, 2, 2.5, 5, 10, ...) com no maximo max_marcas marcacoes"""
    for m in (1, 2, 2.5, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000):
        if faixa / (passo * m) + 1 <= max_marcas + 1e-9:
            return passo * m
    return passo * 1000


def amplia(lim, dados, margem=0.0):
    """lim = (min, max, passo) do Ceron; dados = (min, max) do que sera desenhado. Devolve (min, max, passo, ampliado):
    a faixa do Ceron se os dados cabem nela; senao a uniao das duas, arredondada para fora no passo (o do Ceron, ou um
    multiplo redondo dele se as marcacoes passarem de 11)."""
    a, b, p = lim
    lo, hi = dados
    if lo is None or hi is None or not (math.isfinite(lo) and math.isfinite(hi)):
        return a, b, p, False
    lo, hi = lo - margem * (hi - lo), hi + margem * (hi - lo)
    if lo >= a - 1e-12 and hi <= b + 1e-12:
        return a, b, p, False
    novo_a, novo_b = min(a, lo), max(b, hi)
    if p is None:  # eixo log: decadas inteiras
        return 10 ** math.floor(math.log10(novo_a)), 10 ** math.ceil(math.log10(novo_b)), None, True
    q = _passo_bom(p, novo_b - novo_a)
    novo_a = math.floor(novo_a / q + 1e-9) * q
    novo_b = math.ceil(novo_b / q - 1e-9) * q
    # a faixa do Ceron continua dentro (o arredondamento so cresce)
    novo_a, novo_b = min(novo_a, a), max(novo_b, b)
    return novo_a, novo_b, q, True
