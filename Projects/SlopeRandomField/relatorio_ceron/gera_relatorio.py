#!/usr/bin/env python3
"""gera_relatorio.py: o relatorio inteiro, das pastas de dados ao PDF, no formato da secao 6 de Vargas Ceron et al.
(2025, IJNAMG 49:2459, doi 10.1002/nag.3993).

uso: python3 -I gera_relatorio.py <raiz_dados> <saida> [--etapas=dados,tabelas,figuras_estat,figuras_campos,documento]

Etapas, nesta ordem:
  dados          dados.calcula(<raiz_dados>) e dados.grava(est, <saida>) -> <saida>/estatisticas.json
  tabelas        tabelas.gera(est, <saida>) -> <saida>/tabelas/*.tex
  figuras_estat  figuras_estat.gera(est, <saida>) -> <saida>/figuras/<solo>_fig25..33 e legendas_estat.json
  figuras_campos figuras_campos.gera(est, <raiz_dados>, <saida>) -> <saida>/figuras/<solo>_fig12..14, 16, 24,
                 legendas_campos.json e modos_ruptura.json
  documento      documento.gera(est, <saida>) -> <saida>/relatorio.tex e relatorio.pdf

Os dados so sao lidos. Um erro em dados (uma fusao de pastas com parametros.txt diferentes, por exemplo) para tudo. Um
erro em tabelas ou numa das etapas de figuras e mostrado, e as etapas seguintes rodam: o documento compila com um aviso
visivel no lugar de cada arquivo ausente; o codigo de saida e 1 se alguma etapa falhou. Sem a etapa dados (--etapas sem
ela), o estatisticas.json que ja esta em <saida> e usado."""

import importlib
import json
import os
import sys
import time
import traceback

_AQUI = os.path.dirname(os.path.abspath(__file__))
_PAI = os.path.dirname(_AQUI)  # Projects/SlopeRandomField, onde esta sensibilidade_mc.py
for _p in (_PAI, _AQUI):
    if _p not in sys.path:
        sys.path.insert(0, _p)

os.environ.setdefault("MPLBACKEND", "Agg")

ETAPAS = ("dados", "tabelas", "figuras_estat", "figuras_campos", "documento")


def _uso(msg=None):
    if msg:
        print(msg, file=sys.stderr)
    print(__doc__, file=sys.stderr)
    return 2


def _etapa(nome, func):
    """roda uma etapa, mostra o tempo; devolve (ok, resultado)"""
    t0 = time.time()
    print(f"== {nome}", flush=True)
    try:
        r = func()
    except Exception:
        traceback.print_exc()
        print(f"== {nome}: ERRO ({time.time() - t0:.1f} s)", file=sys.stderr, flush=True)
        return False, None
    print(f"== {nome}: ok ({time.time() - t0:.1f} s)", flush=True)
    return True, r


def _mostra(r):
    if isinstance(r, dict):
        for f in r.get("figuras", []) or []:
            print("   " + str(f))
        for a in r.get("avisos", []) or []:
            print("   AVISO: " + str(a))
    elif isinstance(r, (list, tuple)):
        for f in r:
            print("   " + str(f))
    elif r:
        print("   " + str(r))


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    opcoes = [a for a in argv if a.startswith("--")]
    if len(args) != 2:
        return _uso()
    etapas = list(ETAPAS)
    for o in opcoes:
        if o.startswith("--etapas="):
            etapas = [e.strip() for e in o.split("=", 1)[1].split(",") if e.strip()]
            desconhecidas = [e for e in etapas if e not in ETAPAS]
            if desconhecidas:
                return _uso("etapas desconhecidas: " + ", ".join(desconhecidas))
        elif o in ("--ajuda", "--help"):
            return _uso()
        else:
            return _uso("opção desconhecida: " + o)
    raiz = os.path.abspath(args[0])
    saida = os.path.abspath(args[1])
    if not os.path.isdir(raiz):
        return _uso(f"raiz dos dados não existe: {raiz}")
    os.makedirs(saida, exist_ok=True)
    falhas = []

    # dados: sem eles nada mais faz sentido; um erro aqui para tudo
    if "dados" in etapas:
        import dados
        t0 = time.time()
        print("== dados", flush=True)
        est = dados.calcula(raiz)
        arq = dados.grava(est, saida)
        print(f"   {arq}\n== dados: ok ({time.time() - t0:.1f} s)", flush=True)
    else:
        arq = os.path.join(saida, "estatisticas.json")
        if not os.path.isfile(arq):
            return _uso(f"sem a etapa dados, é preciso {arq}")
        with open(arq, encoding="utf-8") as f:
            est = json.load(f)
        print(f"== dados: lido {arq}", flush=True)

    for nome in ("tabelas", "figuras_estat", "figuras_campos"):
        if nome not in etapas:
            continue
        try:
            mod = importlib.import_module(nome)
        except Exception:
            traceback.print_exc()
            print(f"== {nome}: ERRO ao importar o módulo", file=sys.stderr, flush=True)
            falhas.append(nome)
            continue
        if nome == "figuras_campos":
            ok, r = _etapa(nome, lambda: mod.gera(est, raiz, saida))
        else:
            ok, r = _etapa(nome, lambda: mod.gera(est, saida))
        _mostra(r)
        if not ok:
            falhas.append(nome)

    if "documento" in etapas:
        import documento
        ok, pdf = _etapa("documento", lambda: documento.gera(est, saida))
        _mostra(pdf)
        if not ok:
            falhas.append("documento")

    if falhas:
        print("etapas com erro: " + ", ".join(falhas), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
