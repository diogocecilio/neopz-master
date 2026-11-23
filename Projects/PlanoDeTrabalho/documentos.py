#!/usr/bin/env python3
"""documentos.py: a lista unica dos documentos do trabalho. Dela saem o indice do README.md desta pasta e o "Mapa dos
documentos" do relatorio de situacao (situacao/gera_situacao.py importa mapa_latex), entao os dois nao se desencontram.

Para mudar o indice: edite DOCUMENTOS e rode
    python3 -I documentos.py
O numero de paginas sai do proprio PDF (pdfinfo); um documento que nao existe aparece como "falta".
"""
import glob
import os
import re
import subprocess

AQUI = os.path.dirname(os.path.abspath(__file__))  # Projects/PlanoDeTrabalho
NEOPZ = os.path.dirname(os.path.dirname(AQUI))
README = os.path.join(AQUI, "README.md")
INICIO, FIM = "<!-- indice: gerado por documentos.py; nao edite a mao -->", "<!-- fim do indice -->"

# (grupo, [(caminho a partir da raiz do neopz-master, ou ~/..., descricao)]); a descricao e texto simples (sem LaTeX)
DOCUMENTOS = [
    ("Formulação: o solo não saturado e o fator de segurança (docs/)", [
        ("docs/nao_saturado/biot_nao_saturado_agua_ar5.pdf",
         "formulação poroelástica não saturada (esqueleto, água e ar), a implementação no NeoPZ, Rahardjo et al. (2007) "
         "e benchmarks; a versão 4 está ao lado"),
        ("docs/nao_saturado/formulacao/formulacao.pdf",
         "formulação acoplada u–pw–pa com Mohr–Coulomb associado: discretização, resíduos e jacobiano"),
        ("docs/nao_saturado/resumo/resumo_nao_saturado.pdf", "resumo da formulação u–pw–pa"),
        ("docs/nao_saturado/relatorio_discussao.pdf",
         "contraexemplos a Rahardjo et al. (2007): discussão dos resultados, argumentos e verificação"),
        ("docs/nao_saturado/artigo/artigo.pdf",
         "artigo: a intensidade crítica de chuva em taludes não saturados, e quanto as análises desacopladas e com ar "
         "passivo a superestimam"),
        ("docs/talude_srm/slope_mohr_coulomb/relatorio.pdf",
         "talude Mohr–Coulomb no NeoPZ: correções da plasticidade, FS e rebaixamento"),
        ("docs/talude_srm/slope_mohr_coulomb/nota_uqp.pdf", "nota: formulação de três campos u–q–p"),
        ("docs/talude_srm/fs_arc_length/relatorio_fs.pdf",
         "FS de talude não saturado sob chuva: SRM e GIM com e sem arc-length, sensibilidade e parâmetros"),
    ]),
    ("Formulação: campos aleatórios, FORM e amostragem por importância (docs/confiabilidade_campos)", [
        ("docs/confiabilidade_campos/formulacao/formulacao.pdf",
         "o modelo, as duas bases de Karhunen–Loève (analítica e Nyström) e a mudança de base, o Monte Carlo, a AI, a "
         "variável de controle, o FORM com o gradiente modal (SRM e GIM) e o line sampling (também com o tempo, para a "
         "chuva); nos apêndices, exemplos em Mathematica (exemplo_1d.wl e exemplo_ls.wl)"),
        ("docs/confiabilidade_campos/implementacao_form_ai/implementacao_form_ai.pdf",
         "como o FORM e a AI foram implementados no NeoPZ"),
        ("docs/confiabilidade_campos/ai_nystrom/ai_nystrom.pdf",
         "AI com a KL de Nyström, um roteiro testado em Mathematica"),
        ("docs/confiabilidade_campos/confiabilidade/confiabilidade.pdf",
         "FORM, Monte Carlo e AI em variáveis: a revisão de um código Mathematica e um talude infinito sob chuva (os "
         "números saem de Projects/Confiabilidade)"),
        ("docs/confiabilidade_campos/kl_form_placas/relatorio.pdf",
         "estudo do Ghanem & Spanos (Seção 5.3): KL analítica, Nyström e Galerkin, Monte Carlo e FORM nas placas "
         "quadrada e curva"),
        ("docs/confiabilidade_campos/guia_estudo/guia_estudo.pdf",
         "guia de estudo que liga tudo (inclui o FORM com a base da KL, Eq. 14)"),
    ]),
    ("Resultados: o plano de trabalho (Projects/PlanoDeTrabalho)", [
        ("Projects/PlanoDeTrabalho/README.md",
         "o plano de trabalho (etapas 0 a 5, 2b e 2c), as regras dos documentos e este índice"),
        ("Projects/PlanoDeTrabalho/situacao/situacao.pdf",
         "a situação: o plano, os resultados principais (com as seções do relatório das etapas), as hipóteses e os "
         "próximos passos"),
        ("Projects/PlanoDeTrabalho/ETAPA0/relatorio_etapa0.pdf",
         "etapa 0: a AI da chuva no fermi e a prévia do Cho2010 (SRM antigo; superado pelo das etapas 1 a 3)"),
        ("Projects/PlanoDeTrabalho/ETAPAS1a3/relatorio_etapas.pdf",
         "todos os resultados das etapas 0 a 3 e a referência da chuva: SRM corrigido, AI contra o direto, teste linear, "
         "as duas bases, o gradiente modal (SRM e GIM), o viés, a variável de controle, a fase 4 parcial, o piloto e o "
         "custo"),
        ("Projects/PlanoDeTrabalho/viga_mathematica/relatorio_viga.pdf",
         "a viga com campo aleatório em Mathematica (código em viga_mathematica/codigo): MC, FORM, SORM e ISMC com M, centros baratos (fórmulas de viga), o ponto ótimo do M do FORM e o SORM como função de amostragem"),
        ("Projects/PlanoDeTrabalho/comparacao_metodos/relatorio_comparacao.pdf",
         "comparação com os dados parciais (20/11/2026): Monte Carlo direto, ISMC e line sampling no Cho (r = -0,5 e "
         "-0,7), no rebaixamento do Ceron e na viga; custo para CoV de 5%, convergência, sementes e direções "
         "(comparacao_metodos/gera_relatorio.sh refaz com os dados novos)"),
        ("Projects/PlanoDeTrabalho/contornos_agua/relatorio_contornos.pdf",
         "os contornos da água: Ceron (2025) contra o PZ 'reservatorio' e 'ceron' e a chuva do Rahardjo (2007) contra o "
         "código; esquemas em escala, as duas tabelas de contornos, o Gamma por contorno e domínio, e os próximos passos "
         "(contornos_agua/gera.sh refaz)"),
        ("Projects/PlanoDeTrabalho/ceron_secao6/relatorio_ceron_secao6.pdf",
         "as Seções 6.1 a 6.5 do Ceron (2025) com a água do artigo, figura a figura (Figs. 25 a 33, Tabelas 3 a 6): "
         "FORM modal, line sampling e MC pareado da campanha lsCeronC, com as colunas da água antiga; o de 13/11 está "
         "em arquivo/2026-11-13_relatorio_ceron_agua_antiga.pdf (ceron_secao6/gera.sh refaz depois do "
         "baixar_ls_ceron.sh)"),
        ("Projects/PlanoDeTrabalho/hu2024/relatorio_hu2024.pdf",
         "o caso A da chuva, Hu et al. (2024), na campanha hu2024 (6 cenários, k_s = 2e-5): o determinístico contra a "
         "média de Hu (malhas, os dois k_s, chuva inclinada e vertical), o FORM modal, o line sampling e o MC direto "
         "contra a fragilidade do artigo, o LS contra o MC, por que as P_f de Hu não se reproduzem (o desvio do FS é "
         "0,38 do dele e passa do limite homogêneo) e a primeira passagem (hu2024/gera.sh refaz)"),
        ("Projects/PlanoDeTrabalho/liu2024/relatorio_liu2024.pdf",
         "o caso B da chuva, Liu et al. (2024), na campanha liu2024 (cenários 2 e 3, permanente com campos gaussianos): "
         "o permanente contra a Fig. 3 e o FS contra o Bishop, o FORM modal, o line sampling e o MC direto contra as "
         "Tabelas 2 e 4, o MC que exclui a verificação da Tabela 2, a média do FS (SRM 3,4 % acima do Bishop) que "
         "explica a Tabela 4 e o efeito do campo de k_s (liu2024/gera.sh refaz com as linhas do LS que existirem)"),
        ("Projects/PlanoDeTrabalho/arquivo/", "cópias datadas (AAAA-MM-DD_nome.pdf), nunca editadas"),
    ]),
    ("Resultados: as campanhas", [
        ("Projects/PlanoDeTrabalho/campanhas/relatorio_chuva.pdf", "Monte Carlo da chuva (Rahardjo); a Seção 7 é a AI"),
        ("Projects/PlanoDeTrabalho/campanhas/relatorio_rebaixamento.pdf",
         "rebaixamento do Ceron: Monte Carlo direto e AI (validação: 15,70% contra 15,87%); água antiga (transiente "
         "0,1 h); para a comparação com o Ceron, ver contornos_agua"),
        ("Projects/PlanoDeTrabalho/campanhas/relatorio_solos.pdf", "os três solos do Rahardjo sob chuva"),
        ("Projects/SlopeRandomFieldRefina/relatorio/relatorio.pdf",
         "refinamento por realização e Karhunen–Loève no talude do Ceron (2025); fica com o programa, porque os scripts "
         "dele e os do KLCilindroMika gravam nesta pasta"),
    ]),
    ("Fora do repositório", [
        ("~/Dropbox/UnsaturatedSoil/referencias/", "artigos de terceiros (Cho 2010, Ceron 2025, a tese do Santos, ...)"),
        ("~/projects/runs/LEIAME.md",
         "os dados brutos das rodadas, uma pasta por campanha: o inventário em ordem cronológica (o que é cada pasta, "
         "água e contorno, resultado principal, quem lê e a situação), as três rodadas do Ceron 6.1 a 6.5 e onde está "
         "cada relatório figura a figura"),
        ("~/projects/runs/_arquivo/",
         "as rodadas ultrapassadas e as de estudos concluídos, uma subpasta por data e tema, cada uma com um LEIAME; "
         "os únicos resultados probabilísticos com ar ativo estão em 2026-11-04_chuva_srm_antigo"),
    ]),
]
NOTA = "Como usar cada programa: o README.md da pasta dele, em Projects/<programa>/ (fica com o código)."


def caminho_real(c):
    return os.path.expanduser(c) if c.startswith("~") else os.path.join(NEOPZ, c)


def paginas(c):
    """as paginas do PDF, '' para o que nao e PDF, 'falta' se o arquivo nao existe"""
    p = caminho_real(c)
    if not os.path.exists(p):
        return "falta"
    if not p.endswith(".pdf"):
        return ""
    try:
        s = subprocess.run(["pdfinfo", p], capture_output=True, text=True, timeout=30).stdout
        m = re.search(r"^Pages:\s+(\d+)", s, re.M)
        return m[1] if m else "?"
    except (OSError, subprocess.TimeoutExpired):
        return "?"


def arquivados():
    return sorted(os.path.basename(a) for a in glob.glob(os.path.join(AQUI, "arquivo", "*.pdf")))


def descricao(c, d):
    if c.rstrip("/").endswith("PlanoDeTrabalho/arquivo") and arquivados():
        d += ": " + ", ".join(arquivados())
    return d


def latex(s):
    """texto simples -> LaTeX (o preambulo do relatorio declara o ≈)"""
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\,\%"), ("_", r"\_"), ("#", r"\#"), ("–", "--"),
                 ("−", "$-$"), ("×", r"$\times$"), ("′", "'")):
        s = s.replace(a, b)
    return s


def mapa_latex(este="Projects/PlanoDeTrabalho/situacao/situacao.pdf"):
    """o "Mapa dos documentos" do relatorio de situacao (longtable); <este> e o proprio relatorio (sem paginas: o PDF
    ainda vai ser compilado)"""
    lin = [r"\begin{longtable}{p{0.44\textwidth}rp{0.46\textwidth}}", r"\toprule", r"Documento & p. & O que tem \\",
           r"\midrule", r"\endhead"]
    for grupo, docs in DOCUMENTOS:
        lin.append(r"\multicolumn{3}{p{0.96\textwidth}}{\textbf{" + latex(grupo) + r"}} \\")
        for c, d in docs:
            p, d = ("---", descricao(c, d) + " (este)") if c == este else (paginas(c), descricao(c, d))
            lin.append(r"\path{" + c + "} & " + p + " & " + latex(d) + r" \\")
    lin += [r"\bottomrule", r"\end{longtable}", latex(NOTA)]
    return "\n".join(lin) + "\n"


def indice_markdown():
    lin = [INICIO, ""]
    for grupo, docs in DOCUMENTOS:
        lin += [f"### {grupo}", "", "| Documento | Páginas | O que tem |", "|---|---:|---|"]
        lin += [f"| `{c}` | {paginas(c)} | {descricao(c, d)} |" for c, d in docs]
        lin.append("")
    lin += [NOTA, "", FIM]
    return "\n".join(lin)


def atualiza_readme():
    t = open(README, encoding="utf-8").read()
    i, f = t.find(INICIO), t.find(FIM)
    if i < 0 or f < i:
        raise SystemExit(f"{README}: faltam os marcadores do indice")
    novo = t[:i] + indice_markdown() + t[f + len(FIM):]
    if novo != t:
        open(README, "w", encoding="utf-8").write(novo)
    return novo != t


if __name__ == "__main__":
    mudou = atualiza_readme()
    faltam = [c for _, docs in DOCUMENTOS for c, _ in docs if paginas(c) == "falta"]
    print(f"indice do README {'atualizado' if mudou else 'sem mudancas'}: "
          f"{sum(len(d) for _, d in DOCUMENTOS)} documentos" + (f"; FALTAM: {', '.join(faltam)}" if faltam else ""))
