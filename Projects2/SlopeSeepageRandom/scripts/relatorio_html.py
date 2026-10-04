#!/usr/bin/env python3
"""Relatório comparativo (HTML) a partir de analise_artigo.py.

Uso: relatorio_html.py <saída de analise_artigo.py> <relatorio.html> [previsao.json]
As figuras PNG da análise são embutidas (data URI); o HTML é autocontido.
"""
import base64
import html
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analise_artigo import ARTIGO_DET, ARTIGO_MC, ARTIGO_PF_CHO, GRUPOS, NAO_CONV  # noqa: E402
from campanha_artigo import CASOS  # noqa: E402

ROTULO = {
    "ref": "referência", "covk0": "CoV(k_v) = 0", "covk75": "CoV(k_v) = 75 %", "covk90": "CoV(k_v) = 90 %",
    "covk100": "CoV(k_v) = 100 %", "covc10": "CoV(c) = 10 %", "covc50": "CoV(c) = 50 %", "covc70": "CoV(c) = 70 %",
    "covphi5": "CoV(φ) = 5 %", "covphi15": "CoV(φ) = 15 %", "covphi20": "CoV(φ) = 20 %", "s1.5": "s = 1,5",
    "s2": "s = 2", "s5": "s = 5", "s10": "s = 10", "s20": "s = 20", "s400": "s = 400", "alfa2": "α = 2",
    "alfa3": "α = 3", "alfa4": "α = 4", "alfa5": "α = 5", "hw0.5": "h_w/H = 0,5", "hw0.6": "h_w/H = 0,6",
    "hw0.7": "h_w/H = 0,7", "hw0.8": "h_w/H = 0,8", "hw0.9": "h_w/H = 0,9", "cho_coesivo": "Cho, coesivo",
    "cho_cphi": "Cho, c-φ",
}


def n(v, d=3, pct=False):
    if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
        return "—"
    if pct:
        v *= 100
    s = f"{v:.{d}f}"
    return s.replace(".", ",")


def img(path, alt):
    if not os.path.exists(path):
        return f'<p class="falta">Figura ainda não gerada ({html.escape(os.path.basename(path))}).</p>'
    b = base64.b64encode(open(path, "rb").read()).decode()
    return f'<figure><img src="data:image/png;base64,{b}" alt="{html.escape(alt)}"><figcaption>{alt}</figcaption></figure>'


def diff(a, b):
    """diferença relativa (a/b - 1) em %, com sinal"""
    if a is None or b in (None, 0):
        return "—"
    d = (a / b - 1) * 100
    return ("+" if d >= 0 else "−") + n(abs(d), 1) + " %"


def main(dirres, saida, extra=None):
    R = json.load(open(os.path.join(dirres, "resultados.json")))
    E = json.load(open(extra)) if extra and os.path.exists(extra) else {}
    P = E.get("previsao")
    mc, det = R["mc"], R["det"]
    conv = det.get("conv", {})
    H = []
    w = H.append

    # ------------------------------------------------------------------ cabeçalho e resumo
    nmin = min((mc[c]["N"] for c in mc), default=0)
    nmax = max((mc[c]["N"] for c in mc), default=0)
    ntot = sum(mc[c]["N"] for c in mc)
    milhar = lambda v: f"{v:,}".replace(",", ".")  # noqa: E731
    faixa = milhar(nmin) if nmin == nmax else f"{milhar(nmin)}–{milhar(nmax)}"
    w(f"""<header class="topo">
<p class="sobre">NeoPZ · Projects2/SlopeSeepageRandom · comparação com Vargas Ceron, Cecílio, Linn &amp; Maghous (IJNAMG 2025)</p>
<h1>Talude sob rebaixamento: artigo × elementos finitos</h1>
<p class="lead">Exemplos do artigo recalculados com o código do NeoPZ (elastoplasticidade Mohr-Coulomb incremental, Darcy
por elementos finitos e campos aleatórios de Karhunen-Loève), comparados com os valores publicados.
Monte Carlo: {milhar(ntot)} amostras ({faixa} por caso).</p>
</header>""")
    if E.get("resumo"):
        w('<section id="resumo"><h2>Resumo</h2><ul class="resumo">')
        for item in E["resumo"]:
            w(f"<li>{item}</li>")
        w("</ul></section>")

    # ------------------------------------------------------------------ previsão
    if P:
        w('<section id="previsao"><h2>Previsão de tempo</h2>')
        w(f'<p>{P["texto"]}</p>')
        w('<div class="tab"><table><thead><tr>' + "".join(
            f'<th class="{"num" if i else ""}">{c}</th>' for i, c in enumerate(P["colunas"])) + '</tr></thead><tbody>')
        for linha in P["linhas"]:
            w("<tr>" + "".join(f'<td class="{"num" if i else ""}">{v}</td>' for i, v in enumerate(linha)) + "</tr>")
        w("</tbody></table></div>")
        if P.get("notas"):
            w("<ul>" + "".join(f"<li>{x}</li>" for x in P["notas"]) + "</ul>")
        w("</section>")

    # ------------------------------------------------------------------ como foi calculado
    w("""<section id="metodo"><h2>O que foi comparado</h2>
<p>O artigo calcula o fator de estabilidade Γ pela análise limite cinemática (mecanismos log-espirais, limite
superior), com as forças de percolação −grad u de uma solução semianalítica (K·v'<sub>opt</sub>) ou de elementos
finitos em poropressão (−grad u'<sub>FE</sub>). O código obtém o mesmo Γ como multiplicador das cargas
γ'e<sub>y</sub> − grad u no colapso de uma análise elastoplástica incremental (Mohr-Coulomb associado, quadriláteros
de ordem 2, colapso = não convergência do Newton com bissecção até 0,5 %), com u de Darcy por elementos finitos
(H1, ordem 2). Para a plasticidade perfeita associada, a carga de colapso é única; o elemento finito em
deslocamentos converge para ela de cima com o refinamento, e o mecanismo log-espiral do artigo é um limite superior.
As comparações válidas são, portanto, com as curvas −grad u'<sub>FE</sub> do artigo.</p>
<dl class="pares">
<dt>Peso específico da água</dt><dd>γ<sub>w</sub> = 9,81 kN/m³. O artigo não informa γ<sub>w</sub>, mas os valores da
Fig. 8 em h<sub>w</sub>/H = 0 só são reproduzidos com 9,81 (com 10 a diferença é de 2,3 a 2,5 %).</dd>
<dt>Malha do Monte Carlo</dt><dd>h = 1 m (H = 5 m) com dois níveis de refinamento guiados pelo mecanismo do
problema médio: 4 914 equações no caso de referência, 18 014 no Cho coesivo e 3 194 no Cho c-φ (h = 2 m, H = 10 m).
Os determinísticos usam até cinco níveis.</dd>
<dt>Campos aleatórios</dt><dd>KL de Galerkin (quadriláteros de 9 nós, h<sub>KL</sub> = 1 m) com todos os modos
da KL discreta, lognormais, covariância exponencial; erro de discretização da variância ε<sub>M</sub> ≈ 3,6 % para
(L<sub>x</sub>, L<sub>y</sub>) = (20, 2) m (0,007–2,3 % nos casos da Tabela 4 do artigo), compensado ponto a ponto (o
artigo usa M = 2000 termos com erro &lt; 6 %, sem compensação). c, φ e k<sub>v</sub> independentes.</dd>
<dt>Amostragem</dt><dd>Monte Carlo direto, como no artigo; Pf = P(Γ &lt; 1). Cada amostra é reproduzível isoladamente
(semente, amostra, campo) e a campanha é retomável.</dd>
</dl></section>""")

    # ------------------------------------------------------------------ determinísticos
    w('<section id="det"><h2>Análises determinísticas</h2>')
    # Cho e referência
    cc = {a: conv.get(f"chocoes_h1_a{a}", {}) for a in range(6)}
    cp = {a: conv.get(f"chocphi_h2_a{a}", {}) for a in range(6)}
    rf = {a: conv.get(f"ref_h1_a{a}", {}) for a in range(6)}
    w('<h3>Convergência com a malha</h3><p>Γ e FS do problema médio para níveis crescentes de refinamento adaptativo '
      '(h = 1 m, ou 2 m no exemplo c-φ com H = 10 m). O FE decresce com o refinamento; o valor do artigo é um limite '
      'superior.</p>')
    w('<div class="tab"><table><thead><tr><th>exemplo</th><th>grandeza</th>'
      + "".join(f'<th class="num">nível {a}</th>' for a in range(6))
      + '<th class="num">artigo</th><th class="num">referência</th></tr></thead><tbody>')
    for nome, dd, gk, art, refv in (("Cho coesivo (2:1, c<sub>u</sub> = 23 kPa)", cc, "gamma", 1.354, "F<sub>s</sub> = 1,356 (Cho)"),
                                    ("Cho coesivo", cc, "fs", 1.354, ""),
                                    ("Cho c-φ (1:1, H = 10 m)", cp, "gamma", 1.777, "Chen (1975): 1,777"),
                                    ("Cho c-φ", cp, "fs", 1.203, "F<sub>s</sub> = 1,204 (Cho)"),
                                    ("referência com percolação (Tab. 2)", rf, "gamma", 1.336, ""),
                                    ("referência", rf, "fs", None, "")):
        w(f'<tr><td>{nome}</td><td>{"Γ" if gk == "gamma" else "FS"}</td>'
          + "".join(f'<td class="num">{n(dd[a].get(gk))}</td>' for a in range(6))
          + f'<td class="num">{n(art)}</td><td>{refv}</td></tr>')
    w("</tbody></table></div>")
    extra = []
    if conv.get("ref_gw10_a3", {}).get("gamma"):
        extra.append(f'γ<sub>w</sub> = 10 em vez de 9,81 (nível 3): Γ = {n(conv["ref_gw10_a3"]["gamma"])} '
                     f'(9,81: {n(rf[3].get("gamma"))})')
    if conv.get("ref_dom_a3", {}).get("gamma"):
        extra.append(f'domínio 105 × 55 m em vez de 25 × 10 m (nível 3): Γ = {n(conv["ref_dom_a3"]["gamma"])}')
    for a in range(4):
        g = conv.get(f"ref_h05_a{a}", {}).get("gamma")
        if g:
            extra.append(f"h = 0,5 m, nível {a}: Γ = {n(g)}")
    if extra:
        w("<p>Caso de referência, outras variações: " + "; ".join(extra) + ".</p>")
    for nota in E.get("notas_conv", []):
        w(f"<p>{nota}</p>")

    # Fig. 8
    w('<h3>Fig. 8 — altura crítica com rebaixamento (Wu et al.)</h3>')
    w(img(os.path.join(dirres, "fig8_hcrit.png"),
          "H_crit = Γ·H (problema autossemelhante). Cinza: curvas do artigo (contínua r_p = 0,25 de Wu et al.; "
          "tracejada K·v'_opt; traço-ponto −grad u'_FE). Azul: FE com os parâmetros que reproduzem o artigo em "
          "h_w/H = 0; laranja: FE com os parâmetros da Tabela 1 para o rótulo do painel."))
    cmp8 = det.get("fig8_cmp", [])
    xs = (0., 0.2, 0.4, 0.6, 0.8, 1.0)
    w('<div class="tab"><table><thead><tr><th>painel do artigo (solo do cálculo)</th><th>β</th><th></th>'
      + "".join(f'<th class="num">h<sub>w</sub>/H = {n(x, 1)}</th>' for x in xs) + '</tr></thead><tbody>')
    for panel, solo, betas in (("London Clay", "B: c = 11,7, φ = 24,7°", (30, 60)),
                               ("Israeli Clay", "A: c = 6, φ = 32°", (35, 60))):
        for b in betas:
            pts = {round(r["hwH"], 2): r for r in cmp8 if r["painel"] == panel and r["beta"] == b}
            linhas = (("FE", lambda r: "n.c." if r["nao_conv"] else n(r["fe"], 1)),
                      ("artigo, −grad u'<sub>FE</sub>", lambda r: n(r["grad_u_fe"], 1)),
                      ("FE/artigo − 1", lambda r: "—" if r["nao_conv"] else diff(r["fe"], r["grad_u_fe"])))
            for k, (rot, fn) in enumerate(linhas):
                cab = f"<td>{panel} ({solo})</td><td>{b}°</td>" if k == 0 else "<td></td><td></td>"
                w(f"<tr>{cab}<td>{rot}</td>" + "".join(
                    f'<td class="num">{fn(pts[x]) if x in pts else "…"}</td>' for x in xs) + "</tr>")
    w("</tbody></table></div>")
    w('<p class="nota">H<sub>crit</sub> em m (γ = 18 kN/m³, γ<sub>w</sub> = 9,81 kN/m³, α = 1, malha com 3 níveis). '
      '“n.c.”: sem convergência com a malha (φ próximo de β, sem rebaixamento). Valores do artigo interpolados nas '
      'curvas extraídas do PDF.</p>')
    for nota in E.get("notas_fig8", []):
        w(f"<p>{nota}</p>")

    # Fig. 9 e Tabelas 5/6 (determinístico)
    w('<h3>Fig. 9 — Γ × inclinação e anisotropia</h3>')
    w(img(os.path.join(dirres, "fig9_gamma_beta.png"),
          "Γ × β para α = 1, 5, 10 (H = h_w = 5 m, c = 10 kPa, φ = 30°, γ = 20 kN/m³). Cinza: artigo "
          "(contínua −grad u'_FE, tracejada K·v'_opt, digitalizadas da figura). Azul: FE, nível 3."))
    cmp9 = [r for r in det.get("fig9_cmp", []) if r["curva"] == "grad"]
    if cmp9:
        bs = (30, 45, 60, 75, 90)
        w('<div class="tab"><table><thead><tr><th>α</th><th></th>' + "".join(f'<th class="num">β = {b}°</th>' for b in bs)
          + '</tr></thead><tbody>')
        for a9 in (1, 5, 10):
            pts = {r["beta"]: r for r in cmp9 if r["alfa"] == a9}
            for k, (rot, fn) in enumerate((("FE", lambda r: n(r["fe"])), ("artigo, −grad u'<sub>FE</sub>", lambda r: n(r["artigo"])),
                                           ("FE/artigo − 1", lambda r: diff(r["fe"], r["artigo"])))):
                w(f'<tr><td>{a9 if k == 0 else ""}</td><td>{rot}</td>' + "".join(
                    f'<td class="num">{fn(pts[b]) if b in pts else "—"}</td>' for b in bs) + "</tr>")
        w("</tbody></table></div>")
    for nota in E.get("notas_fig9", []):
        w(f"<p>{nota}</p>")
    t56 = det.get("tab56", {})
    w('<h3>Tabelas 5 e 6 — coluna determinística</h3><div class="tab"><table><thead><tr><th>caso</th>'
      '<th class="num">Γ FE (nível 2, malha do MC)</th><th class="num">Γ FE (nível 3)</th><th class="num">Γ artigo</th>'
      '<th class="num">FE nível 3 / artigo</th></tr></thead><tbody>')
    linhas = [("α = 1 (referência)", "ref", conv.get("ref_h1_a2", {}), conv.get("ref_h1_a3", {}))]
    for a in (2, 3, 4, 5):
        linhas.append((f"α = {a}", f"alfa{a}", t56.get(f"alfa{a}_a2", {}), t56.get(f"alfa{a}", {})))
    for h, hw in ((0.5, "2.5"), (0.6, "3"), (0.7, "3.5"), (0.8, "4"), (0.9, "4.5")):
        linhas.append((f"h_w/H = {n(h, 1)}", f"hw{h}", t56.get(f"hw{hw}_a2", {}), t56.get(f"hw{hw}", {})))
    for nome, c, d2, d3 in linhas:
        a = ARTIGO_DET.get(c)
        w(f'<tr><td>{nome}</td><td class="num">{n(d2.get("gamma"))}</td><td class="num">{n(d3.get("gamma"))}</td>'
          f'<td class="num">{n(a)}</td><td class="num">{diff(d3.get("gamma"), a)}</td></tr>')
    w("</tbody></table></div>")

    # Fig. 5
    w('<h3>Fig. 5 — funcional hidráulico</h3>')
    w(img(os.path.join(dirres, "fig5_funcional.png"),
          "J(u)/(k_h H² γw²) do FE em poropressão (malha convergida) para dois domínios, entre as estimativas do "
          "artigo J(u'_FE) (superior) e −J*(v'_opt) (inferior). O valor depende do tamanho do domínio, que o artigo "
          "não informa para o FE; com o domínio de 10H do modelo semianalítico (L_m = 10H) ele fica entre os dois."))
    w("</section>")

    # ------------------------------------------------------------------ Monte Carlo
    w('<section id="mc"><h2>Monte Carlo</h2>')
    w('<p>μ, σ e CoV de Γ; Pf com intervalo de confiança de 95 % (Wilson). As colunas com * reescalam cada amostra por '
      'Γ<sub>artigo</sub>/Γ<sub>FE</sub> do problema médio na mesma malha, isto é, retiram a diferença determinística '
      '(discretização e forças de percolação) e deixam só o efeito da variabilidade.</p>')
    for titulo, casos in GRUPOS:
        w(f"<h3>{html.escape(titulo)}</h3>")
        w('<div class="tab"><table><thead><tr><th>caso</th><th class="num">N</th><th class="num">μ</th>'
          '<th class="num">σ</th><th class="num">CoV</th><th class="num">Pf (IC 95 %)</th><th class="num">CoV(Pf)</th>'
          '<th class="num">μ*</th><th class="num">σ*</th><th class="num">Pf*</th>'
          '<th class="num art">μ artigo</th><th class="num art">σ</th><th class="num art">CoV</th>'
          '<th class="num art">Pf</th></tr></thead><tbody>')
        for c in casos:
            st = mc.get(c)
            if not st:
                continue
            a = st.get("artigo", {})
            e = st.get("escalado") or {}
            pfart = n(a.get("pf"), 2, True) + " %" if "pf" in a else "—"
            if "pf_cho" in a:
                pfart += f' <span class="sub">(Cho {n(a["pf_cho"], 2, True)} %)</span>'
            w(f'<tr><td>{ROTULO.get(c, c)}</td><td class="num">{st["N"]:,}</td><td class="num">{n(st["mu"])}</td>'
              f'<td class="num">{n(st["sd"])}</td><td class="num">{n(st["cov"], 1, True)} %</td>'
              f'<td class="num">{n(st["pf"], 2, True)} % <span class="sub">({n(st["pf_lo"], 1, True)}–'
              f'{n(st["pf_hi"], 1, True)})</span></td><td class="num">{n(st["covpf"], 1, True)} %</td>'
              f'<td class="num">{n(e.get("mu"))}</td><td class="num">{n(e.get("sd"))}</td>'
              f'<td class="num">{n(e.get("pf"), 2, True) + " %" if e else "—"}</td>'
              f'<td class="num art">{n(a.get("mu"))}{"<sup>†</sup>" if "fonte_momentos" in a else ""}</td>'
              f'<td class="num art">{n(a.get("sd"))}</td>'
              f'<td class="num art">{n(a.get("cov"), 1, True) + " %" if "cov" in a else "—"}</td>'
              f'<td class="num art">{pfart}</td></tr>'.replace(f'{st["N"]:,}', f'{st["N"]:,}'.replace(",", ".")))
        w("</tbody></table></div>")
    for nota in E.get("notas_mc", []):
        w(f"<p>{nota}</p>")
    w(img(os.path.join(dirres, "tendencias.png"),
          "Tendências das Tabelas 3–6: Pf, μ(Γ) e CoV(Γ). Vermelho: artigo; azul: FE (barras: IC 95 % de Pf); "
          "verde: FE reescalado pelo Γ determinístico."))
    w(img(os.path.join(dirres, "densidades.png"),
          "Densidade de Γ: referência (lognormal ajustada no artigo, μ = 1,353, σ = 0,318) e exemplos de Cho "
          "(curvas do artigo extraídas das Figs. 17 e 21)."))
    w(img(os.path.join(dirres, "convergencia_pf.png"), "Convergência de Pf com o número de amostras."))
    # Spearman e modos de ruptura
    sp = mc.get("ref", {}).get("spearman")
    if sp and "mecanismo" in sp:
        w('<h3>Sensibilidade — correlação parcial de Spearman de 2ª ordem (Fig. 26)</h3>')
        w('<div class="tab"><table><thead><tr><th>variável</th><th class="num">FE, médias no mecanismo</th>'
          '<th class="num">FE, médias no domínio</th><th class="num">artigo</th></tr></thead><tbody>')
        for k, lab, art in (("c", "c", 0.879), ("phi", "φ", 0.508), ("kv", "k<sub>v</sub>", -0.036)):
            w(f'<tr><td>{lab}</td><td class="num">{n(sp["mecanismo"][k][0])}</td>'
              f'<td class="num">{n(sp["dominio"][k][0])}</td><td class="num">{n(art)}</td></tr>')
        w(f'</tbody></table></div><p class="nota">N = {sp["N"]}. Mecanismo: c e φ médios na banda de cisalhamento '
          '(‖Δε<sup>p</sup>‖ ≥ 10 % do máximo no colapso, ponderados), k<sub>v</sub> médio na massa que se move '
          '(‖Δu‖ ≥ 30 % do máximo), como as médias ao longo da superfície de ruptura e no volume deslizante do artigo.</p>')
    modos = R.get("modos")
    if modos:
        w('<h3>Modos de ruptura (Figs. 16 e 20)</h3><div class="tab"><table><thead><tr><th>exemplo</th>'
          '<th class="num">N</th><th class="num">acima do pé</th><th class="num">pelo pé</th>'
          '<th class="num">abaixo do pé</th><th class="num art">artigo (acima / pelo / abaixo)</th></tr></thead><tbody>')
        art = {"cho_coesivo": "0,55 / 6,76 / 92,68 %", "cho_cphi": "19,20 / 80,80 / 0 %"}
        for c, m in modos.items():
            tot = sum(m.values())
            w(f'<tr><td>{ROTULO.get(c, c)}</td><td class="num">{tot}</td>'
              + "".join(f'<td class="num">{n(m.get(k, 0) / tot, 1, True)} %</td>' for k in ("acima", "pe", "abaixo"))
              + f'<td class="num art">{art.get(c, "—")}</td></tr>')
        w('</tbody></table></div><p class="nota">Classificação pela banda ‖Δε<sup>p</sup>‖ ≥ 20 % do máximo no colapso: '
          'abaixo do pé se desce mais de 0,1H abaixo do nível do pé; pelo pé se passa a menos de 0,1H do pé.</p>')
    w("</section>")

    # ------------------------------------------------------------------ achados e reprodução
    if R.get("achados"):
        w('<section id="achados"><h2>Observações sobre o artigo</h2><ol>')
        for a in R["achados"]:
            w(f"<li>{a}</li>")
        w("</ol></section>")
    w("""<section id="reproduzir"><h2>Como reproduzir</h2>
<p>Branch <code>claude/great-clarke-xist30</code> do neopz-master; dados, logs e este relatório em
<code>Projects2/SlopeSeepageRandom/resultados/artigo2025/</code>.</p>
<pre><code>sudo apt install cmake g++ make liblapack-dev libblas-dev liblapacke-dev
# compila, gera a fila (blocos de 50 amostras intercalados entre os 28 casos) e roda em segundo plano
bash Projects2/SlopeSeepageRandom/scripts/rodar_campanha.sh 1000      # [alvo: N | cov5 | artigo] [processos] [dir]
S=Projects2/SlopeSeepageRandom/scripts
PROCESSOS=8 python3 $S/campanha_artigo.py status ~/campanha_artigo     # andamento e previsão
python3 $S/campanha_artigo.py pacote ~/campanha_artigo resultados.tar.gz
# análise e relatório (det: logs do comando det; ref: dados extraídos do artigo)
python3 $S/analise_artigo.py campanha det ref saida
python3 $S/relatorio_html.py saida relatorio.html saida/extra.json</code></pre></section>""")

    css = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "relatorio.css")).read()
    pagina = ("<title>Talude com percolação no NeoPZ</title>\n"
              '<link rel="preconnect" href="https://fonts.googleapis.com">'
              '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600'
              '&family=IBM+Plex+Mono:wght@400;500&family=Newsreader:opsz,wght@6..72,500;6..72,600&display=swap">'
              f"<style>{css}</style>\n<main>" + "\n".join(H) + "</main>")
    open(saida, "w").write(pagina)
    print(f"{saida}: {len(pagina) / 1e6:.2f} MB")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    main(*sys.argv[1:])
