#!/usr/bin/env python3
"""ismc_parada.py: a etapa 3b do plano de trabalho (campanha etapa3b do cluster). Quantas realizacoes a ISMC (amostragem
por importancia com o deslocamento no ponto de projeto de um FORM) precisa para o CoV alvo da P_f, por centro (o FORM) e
por r, contra os diretos de referencia.

Para cada caso de ISMC do arquivo de casos (os com cov_alvo=), as realizacoes na ordem (real = 0, 1, ...; so as com
situacao ok e FS_min > 0, como o amostragem_importancia.py) dao a P_f acumulada, P_k = (1/k) sum_{i<=k} I_i w_i, e o seu
CoV, sqrt(var(I w) / k) / P_k. O N do alvo e o primeiro k com CoV <= alvo, k >= cov_n e pelo menos cov_falhas falhas
sorteadas (a regra do tarefas.sh); o N estavel e o k a partir do qual o CoV nao volta a passar do alvo. A P_f no N do
alvo (a que a parada entrega) e a P_f com todas as realizacoes (os blocos que ja rodavam terminaram) vao juntas: a
diferenca mede o vies da parada opcional, que deve ficar dentro do ruido.

Por centro e por r: a media das sementes (no N do alvo), o desvio entre elas, a P_f das sementes juntas e o z contra o
direto, o N do alvo (media, minimo e maximo) e o custo total (as avaliacoes do FORM + o N do alvo). O ganho e o N que o
direto precisaria para o mesmo CoV, (1 - P) / (P alvo^2) com a P do direto, dividido pelo custo.

Uso: python3 -I ismc_parada.py <pasta da campanha> <casos_etapa3b.txt> [r05=<pasta do direto>] [r07=<pasta do direto>]
     [figuras=<pasta>] [conf_ref=<pasta do r07_ai_form20 da etapa3>] [cov_n=<minimo>] [forms=<pasta>]
cov_n= troca o cov_n dos casos (o minimo de realizacoes da regra de parada; o etapa3b_gim roda com cov_n=2000 para
parear, e cov_n=500 da o N alvo comparavel ao da etapa 3b). forms= e a pasta com os form_*/form.json (padrao: a da
campanha; no etapa3b_gim os FORM sao os da etapa 3b).
Os diretos padrao sao os baixados das etapas 1 e 3: ~/projects/runs/fermi_etapa1/etapa1/r05_direto e
~/projects/runs/fermi_etapa3/etapa3/r07_direto."""
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "SlopeRandomFieldRefinement",
                                "posproc"))
import amostragem_importancia as ai  # noqa: E402  (a leitura das partes mc_<a>_<b>.csv)


def serie(pasta):
    """os I w das realizacoes na ordem, ate o primeiro buraco (um bloco que ainda nao acabou)"""
    rs = ai.linhas(pasta)
    excl = ai.EXCLUIDAS.get(pasta, {})
    out, n = [], 0
    while n in rs or n in excl:
        if n in rs:
            r = rs[n]
            out.append(math.exp(float(r.get("log_peso") or 0.)) * (float(r["FS_min"]) < 1.))
        n += 1
    return out, len(excl)


def acumulada(iw):
    """P_k e CoV_k para k = 1..n"""
    s1 = s2 = 0.
    pk, cvk = [], []
    for k, x in enumerate(iw, 1):
        s1 += x
        s2 += x * x
        p = s1 / k
        if k > 1 and p > 0:
            v = max((s2 - k * p * p) / (k - 1), 0.)
            cvk.append(math.sqrt(v / k) / p)
        else:
            cvk.append(float("inf"))
        pk.append(p)
    return pk, cvk


def estatisticas(iw):
    n = len(iw)
    s1 = sum(iw)
    s2 = sum(x * x for x in iw)
    p = s1 / n if n else 0.
    cv = math.sqrt(max((s2 - n * p * p) / (n - 1), 0.) / n) / p if n > 1 and p > 0 else float("inf")
    ess = s1 * s1 / s2 if s2 > 0 else 0.
    maior = max(iw) / s1 if s1 > 0 else float("nan")
    return p, cv, sum(1 for x in iw if x > 0), ess, maior


def direto(pasta):
    rs = ai.linhas(pasta)
    n = len(rs)
    f = sum(1 for r in rs.values() if float(r["FS_min"]) < 1.)
    p = f / n
    return {"N": n, "falhas": f, "Pf": p, "CoV": math.sqrt((1 - p) / (n * p)) if p > 0 else float("inf")}


def le_casos(arq):
    casos = []
    for linha in open(arq):
        linha = linha.strip()
        if not linha or linha.startswith("#"):
            continue
        id_, fase, n, *opc = linha.split()
        o = dict(w.split("=", 1) for w in opc if "=" in w)
        casos.append((id_, int(n), o))
    return casos


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    camp, arq_casos = sys.argv[1], sys.argv[2]
    p = {"r05": os.path.expanduser("~/projects/runs/fermi_etapa1/etapa1/r05_direto"),
         "r07": os.path.expanduser("~/projects/runs/fermi_etapa3/etapa3/r07_direto"), "figuras": ""}
    for a in sys.argv[3:]:
        k, v = a.split("=", 1)
        p[k] = v
    refs = {}
    for r in ("r05", "r07"):
        if os.path.isdir(p[r]):
            refs[r] = direto(p[r])
            d = refs[r]
            print(f"direto {r}: N = {d['N']}, {d['falhas']} falhas, P_f = {100 * d['Pf']:.4f}% (CoV {100 * d['CoV']:.1f}%)"
                  f"  [{p[r]}]")
        else:
            print(f"direto {r}: sem a pasta {p[r]}")
    # os FORM
    forms = {}
    for nome in sorted(os.listdir(p.get("forms", camp))):
        f = os.path.join(p.get("forms", camp), nome, "form.json")
        if nome.startswith("form_") and os.path.exists(f):
            d = json.load(open(f))
            forms[nome[5:]] = d
    if forms:
        print("\nFORM (base de Nystrom, M termos por campo):")
        print(f"  {'FORM':<10} {'beta':>7} {'P_f FORM (%)':>12} {'iter.':>5} {'analises':>8}  convergiu")
        for nome, d in forms.items():
            print(f"  {nome:<10} {d['beta']:7.4f} {100 * d['P_f']:12.4f} {d['iteracoes']:5d} {d['avaliacoes']:8d}  "
                  f"{d['convergiu']}")
        for r in ("r05", "r07"):
            us = {c: forms[f"{r}_{c}"]["u"] for c in ("M20", "M100", "gim20") if f"{r}_{c}" in forms}
            if "M100" in us:
                def emM100(u):  # u (2M: c' e phi') na base de 100 termos por campo
                    m = len(u) // 2
                    return u[:m] + [0.] * (100 - m) + u[m:] + [0.] * (100 - m)
                ref = emM100(us["M100"])
                for c, u in us.items():
                    if c == "M100":
                        continue
                    v = emM100(u)
                    cos = sum(a * b for a, b in zip(v, ref)) / math.sqrt(sum(a * a for a in v) * sum(b * b for b in ref))
                    print(f"  {r}: angulo entre o u* de {c} e o de M100 = {math.degrees(math.acos(min(1., cos))):.1f} graus")
    # as ISMC
    grupos = {}
    print("\nISMC por semente (N alvo: o primeiro k com CoV <= alvo; estavel: daqui em diante sempre <= alvo):")
    print(f"  {'caso':<20} {'N alvo':>6} {'estavel':>7} {'P_f alvo (%)':>12} {'N fim':>6} {'P_f fim (%)':>11} "
          f"{'CoV fim':>7} {'falhas':>6} {'ESS_F':>6} {'maior':>6} {'z fim':>6}  ATINGIU")
    curvas = {}
    for id_, nteto, o in le_casos(arq_casos):
        if "cov_alvo" not in o:
            continue
        m = re.match(r"(r0[57])_nys_(\w+)_s(\d+)$", id_)
        if not m:
            continue
        r, centro, sem = m.group(1), m.group(2), int(m.group(3))
        pasta = os.path.join(camp, id_)
        if not os.path.isdir(pasta):
            continue
        alvo, nmin, fmin = float(o["cov_alvo"]), int(p.get("cov_n", o.get("cov_n", 500))), int(o.get("cov_falhas", 100))
        iw, nexcl = serie(pasta)
        if len(iw) < 2:
            continue
        pk, cvk = acumulada(iw)
        falhas_k, nf = [], 0
        for x in iw:
            nf += x > 0
            falhas_k.append(nf)
        nalvo = next((k for k in range(1, len(iw) + 1) if k >= nmin and falhas_k[k - 1] >= fmin and cvk[k - 1] <= alvo),
                     None)
        nest = None
        if cvk[-1] <= alvo:
            nest = len(iw)
            while nest > 1 and cvk[nest - 2] <= alvo:
                nest -= 1
        pf, cv, f, ess, maior = estatisticas(iw)
        ref = refs.get(r)
        z = (pf - ref["Pf"]) / math.sqrt((cv * pf) ** 2 + (ref["CoV"] * ref["Pf"]) ** 2) if ref else float("nan")
        atingiu = open(os.path.join(pasta, "ATINGIU")).read().split()[:4] if os.path.exists(
            os.path.join(pasta, "ATINGIU")) else "-"
        print(f"  {id_:<20} {nalvo or '-':>6} {nest or '-':>7} "
              f"{(100 * pk[nalvo - 1]) if nalvo else float('nan'):12.4f} {len(iw):6d} {100 * pf:11.4f} {100 * cv:6.1f}% "
              f"{f:6d} {ess:6.0f} {maior:6.3f} {z:+6.2f}  {' '.join(atingiu) if atingiu != '-' else '-'}"
              f"{f'  ({nexcl} excluidas)' if nexcl else ''}")
        grupos.setdefault((r, centro), []).append({"sem": sem, "nalvo": nalvo, "palvo": pk[nalvo - 1] if nalvo else None,
                                                    "iw": iw, "pf": pf, "cv": cv})
        curvas[id_] = (pk, cvk)
    # por centro e por r
    print("\nPor centro e por r (as sementes):")
    print(f"  {'r':<4} {'centro':<6} {'sem.':>4} {'N alvo: media (min-max)':>24} {'P_f alvo: media (dp) %':>24} "
          f"{'P_f juntas (%)':>16} {'z':>6} {'FORM':>5} {'custo':>6} {'ganho':>6}")
    for (r, centro), gs in sorted(grupos.items()):
        ns = [g["nalvo"] for g in gs if g["nalvo"]]
        ps = [g["palvo"] for g in gs if g["palvo"] is not None]
        todos = [x for g in gs for x in g["iw"]]
        pj, cvj, _, _, _ = estatisticas(todos)
        ref = refs.get(r)
        z = (pj - ref["Pf"]) / math.sqrt((cvj * pj) ** 2 + (ref["CoV"] * ref["Pf"]) ** 2) if ref else float("nan")
        nform = forms.get(f"{r}_{centro}", {}).get("avaliacoes", 0)
        if ns:
            nm = sum(ns) / len(ns)
            custo = nform + nm
            alvo = 0.05
            ganho = (1 - ref["Pf"]) / (ref["Pf"] * alvo ** 2) / custo if ref else float("nan")
            pm = sum(ps) / len(ps)
            dp = math.sqrt(sum((x - pm) ** 2 for x in ps) / (len(ps) - 1)) if len(ps) > 1 else float("nan")
            print(f"  {r:<4} {centro:<6} {len(gs):4d} {nm:10.0f} ({min(ns)}-{max(ns)}){'':>7} "
                  f"{100 * pm:11.4f} ({100 * dp:.4f}){'':>3} {100 * pj:8.4f} ({100 * cvj:.1f}%) {z:+6.2f} {nform:5d} "
                  f"{custo:6.0f} {ganho:5.1f}x")
        else:
            print(f"  {r:<4} {centro:<6} {len(gs):4d}   (nenhuma semente chegou ao alvo ainda)  P_f juntas "
                  f"{100 * pj:.4f}% ({100 * cvj:.1f}%), z {z:+.2f}")
    # a conferencia de versao: as mesmas realizacoes (semente 3, mesmo deslocamento) da etapa 3 tem de dar o mesmo FS
    conf = os.path.join(camp, "conf_r07_s3")
    antigo = p.get("conf_ref", os.path.expanduser("~/projects/runs/fermi_etapa3/etapa3/r07_ai_form20"))
    if os.path.isdir(conf) and os.path.isdir(antigo):
        a, b = ai.linhas(conf), ai.linhas(antigo)
        comuns = sorted(set(a) & set(b))
        if comuns:
            dif = max(abs(float(a[n]["FS_min"]) - float(b[n]["FS_min"])) for n in comuns)
            print(f"\nconferencia de versao (conf_r07_s3 contra {antigo}): {len(comuns)} realizacoes, max |dFS| = {dif:.2e}"
                  f" ({'ok' if dif < 1e-3 else 'DIFERENTE: o programa mudou'})")
    if p["figuras"]:
        figuras(p["figuras"], grupos, curvas, refs)


def figuras(pasta, grupos, curvas, refs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(pasta, exist_ok=True)
    for r in ("r05", "r07"):
        cs = sorted(c for (rr, c) in grupos if rr == r)
        if not cs:
            continue
        fig, ax = plt.subplots(2, len(cs), figsize=(4.2 * len(cs), 6.5), squeeze=False)
        for j, c in enumerate(cs):
            for g in grupos[(r, c)]:
                pk, cvk = curvas[f"{r}_nys_{c}_s{g['sem']}"]
                k = range(1, len(pk) + 1)
                ax[0, j].plot(k, [100 * x for x in pk], lw=0.8, label=f"s{g['sem']}")
                ax[1, j].plot(k, [100 * x for x in cvk], lw=0.8)
            if r in refs:
                d = refs[r]
                for v, ls in ((d["Pf"], "-"), (d["Pf"] * (1 - 1.96 * d["CoV"]), ":"), (d["Pf"] * (1 + 1.96 * d["CoV"]), ":")):
                    ax[0, j].axhline(100 * v, color="k", ls=ls, lw=0.8)
            ax[1, j].axhline(5, color="k", lw=0.8)
            ax[0, j].set_title(f"{r}, centro {c}")
            ax[0, j].set_xscale("log"); ax[1, j].set_xscale("log"); ax[1, j].set_yscale("log")
            ax[0, j].set_ylabel("P_f (%)"); ax[1, j].set_ylabel("CoV (%)"); ax[1, j].set_xlabel("realizacoes")
            ax[0, j].legend(fontsize=7)
            ax[1, j].set_ylim(1, 100)
        fig.tight_layout()
        fig.savefig(os.path.join(pasta, f"ismc_{r}.pdf"))
        plt.close(fig)
        print(f"figura: {os.path.join(pasta, f'ismc_{r}.pdf')}")


if __name__ == "__main__":
    main()
