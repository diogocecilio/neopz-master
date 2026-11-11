"""estilo.py: estilo comum das figuras de Karhunen-Loeve do relatorio (paleta categorica de
../../SlopeRandomField/comparacao_ceron/figuras_sat.py, decimal com virgula) e os caminhos do relatorio."""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
PROJETOS = os.path.normpath(os.path.join(AQUI, "..", ".."))
RELATORIO = os.environ.get("RELATORIO", os.path.join(PROJETOS, "SlopeRandomFieldRefina", "relatorio"))
FIGURAS = os.path.join(RELATORIO, "figuras")
TABELAS = os.path.join(RELATORIO, "tabelas")
DADOS = os.path.join(RELATORIO, "dados")
MALHAS = os.path.join(PROJETOS, "LearningMohrCoulomb2", "cluster", "malhas")
SCRATCH = os.environ.get("KL_SCRATCH", "/tmp/claude-0/-home-user-neopz-master/674eb771-9e25-5098-bcbd-eda87ce4d4e9/"
                                       "scratchpad")
for d in (FIGURAS, TABELAS, DADOS):
    os.makedirs(d, exist_ok=True)

PALETA = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
TINTA, TINTA2, GRADE = "#0b0b0b", "#52514e", "#dddcd8"
NEUTRO = "#8a8984"
plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7,
    "axes.titlesize": 8.5, "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2,
    "ytick.color": TINTA2, "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.5, "lines.linewidth": 1.4,
    "legend.frameon": False, "savefig.dpi": 200, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
    "axes.unicode_minus": True, "figure.dpi": 100,
})
LARGURA = 6.5  # polegadas: \linewidth do A4 com margens de 2,2 cm (16,6 cm)

# divergente com o meio cinza neutro (azul da paleta <-> vermelho da paleta)
DIVERGENTE = LinearSegmentedColormap.from_list("azul_cinza_vermelho", ["#1d4f91", "#2a78d6", "#9cbde6", "#ecebe8",
                                                                      "#f1a9a3", "#e34948", "#9e2520"])


def num(v, casas=None, sig=None):
    """numero com virgula decimal"""
    if sig is not None:
        s = "%.*g" % (sig, v)
    elif casas is not None:
        s = "%.*f" % (casas, v)
    else:
        s = "%g" % v
    return s.replace(".", ",").replace("-", "−")


def cientifico(v, sig=2):
    """2,9e-3 -> 2,9e−3 (como no README); para o LaTeX use cientifico_tex"""
    if v == 0 or not np.isfinite(v):
        return "0" if v == 0 else "-"
    m, e = ("%.*e" % (sig - 1, v)).split("e")
    return "%se%d" % (m.replace(".", ","), int(e))


def cientifico_tex(v, sig=2):
    """2,9 \\cdot 10^{-3} em modo matematico do LaTeX"""
    if v == 0:
        return "$0$"
    if not np.isfinite(v):
        return "--"
    m, e = ("%.*e" % (sig - 1, v)).split("e")
    return r"$%s\cdot10^{%d}$" % (m.replace(".", "{,}"), int(e))


def fmt_virgula(casas=None):
    """FuncFormatter com virgula decimal (casas fixas ou %g)"""
    def f(v, _pos):
        if casas is None:
            s = ("%g" % v)
        else:
            s = "%.*f" % (casas, v)
        return s.replace(".", ",").replace("-", "−")
    return FuncFormatter(f)


def fmt_log10():
    """rotulos de eixo log: 10^k como potencia; 1, 10, 100, 0,1, 0,01 por extenso so ate 1000 e ate 0,001"""
    def f(v, _pos):
        if v <= 0:
            return ""
        e = np.log10(v)
        if abs(e - round(e)) > 1e-9:
            return ""
        e = int(round(e))
        if -3 <= e <= 3:
            return ("%g" % v).replace(".", ",")
        return r"$10^{%d}$" % e
    return FuncFormatter(f)


def eixo_log(ax, eixo="both"):
    """eixos logaritmicos com rotulos 10^k e virgula"""
    for nome in ("x", "y") if eixo == "both" else (eixo,):
        a = ax.xaxis if nome == "x" else ax.yaxis
        a.set_major_locator(LogLocator(base=10, numticks=12))
        a.set_major_formatter(fmt_log10())
        a.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=12))
        a.set_minor_formatter(NullFormatter())


def virgula(ax, x=None, y=None):
    """virgula decimal nos dois eixos lineares; x/y: numero de casas (None: %g)"""
    ax.xaxis.set_major_formatter(fmt_virgula(x))
    ax.yaxis.set_major_formatter(fmt_virgula(y))


def contorno_talude(ax, **kw):
    """o contorno do talude (0,0) (25,0) (25,5) (15,5) (10,10) (0,10)"""
    p = np.array([[0., 0.], [25., 0.], [25., 5.], [15., 5.], [10., 10.], [0., 10.], [0., 0.]])
    opt = dict(color=TINTA, lw=0.8)
    opt.update(kw)
    ax.plot(p[:, 0], p[:, 1], **opt)


def salva(fig, nome, png=False, dpi=200):
    """grava em FIGURAS: PDF vetorial (linhas) ou PNG (mapas densos); devolve o caminho"""
    caminho = os.path.join(FIGURAS, nome + (".png" if png else ".pdf"))
    if png:
        fig.savefig(caminho, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
    else:
        fig.savefig(caminho, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    tam = os.path.getsize(caminho) / 1e6
    print("  %s (%.2f MB)" % (caminho, tam))
    if tam > 1.5:
        print("  AVISO: %s passa de 1,5 MB" % caminho)
    return caminho
