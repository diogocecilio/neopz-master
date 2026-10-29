"""Malha estruturada de triângulos de uma coluna (ensaio de drenagem de Liakopoulos 1965), sem o gmsh.

Retângulo [0, largura] x [0, altura] com nx x ny células, cada uma dividida em dois triângulos. Grupos físicos das
linhas, como nas malhas dos taludes: 1 base, 2 lateral esquerda, 3 lateral direita, 5 topo (a superfície do
programa SlopeRainfall). Saída: .msh 2.2 ASCII.

Uso: python3 malha_coluna.py [altura = 1] [largura = 0,1] [nx = 2] [ny = 80] [saida = coluna.msh]
"""
import sys


def gerar(arq, altura=1.0, largura=0.1, nx=2, ny=80):
    no = lambda i, j: j * (nx + 1) + i + 1
    nos = [(no(i, j), largura * i / nx, altura * j / ny) for j in range(ny + 1) for i in range(nx + 1)]
    linhas = []
    for i in range(nx):
        linhas.append((1, no(i, 0), no(i + 1, 0)))              # base
        linhas.append((5, no(i + 1, ny), no(i, ny)))            # topo
    for j in range(ny):
        linhas.append((2, no(0, j + 1), no(0, j)))              # lateral esquerda
        linhas.append((3, no(nx, j), no(nx, j + 1)))            # lateral direita
    tri = []
    for j in range(ny):
        for i in range(nx):
            a, b, c, d = no(i, j), no(i + 1, j), no(i + 1, j + 1), no(i, j + 1)
            tri += [(a, b, c), (a, c, d)]
    with open(arq, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n" % len(nos))
        for k, x, y in nos:
            f.write(f"{k} {x:.10g} {y:.10g} 0\n")
        f.write("$EndNodes\n$Elements\n%d\n" % (len(linhas) + len(tri)))
        k = 0
        for tag, a, b in linhas:
            k += 1
            f.write(f"{k} 1 2 {tag} {tag} {a} {b}\n")
        for a, b, c in tri:
            k += 1
            f.write(f"{k} 2 2 1 1 {a} {b} {c}\n")
        f.write("$EndElements\n")
    return len(tri)


if __name__ == "__main__":
    v = sys.argv[1:]
    altura = float(v[0]) if len(v) > 0 else 1.0
    largura = float(v[1]) if len(v) > 1 else 0.1
    nx = int(v[2]) if len(v) > 2 else 2
    ny = int(v[3]) if len(v) > 3 else 80
    arq = v[4] if len(v) > 4 else "coluna.msh"
    print(arq, gerar(arq, altura, largura, nx, ny), "triângulos")
