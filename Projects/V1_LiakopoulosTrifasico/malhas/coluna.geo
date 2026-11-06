// Coluna de areia de Liakopoulos (Secao 10 do resumo "Solos nao saturados"; Callari e Abati 2009, Fig. 2):
// retangulo de 0,10 m x 1,00 m com 4 x 40 quadrados de 0,025 m, cada um dividido em dois triangulos
// (320 elementos). Unidades: m. y para cima, base em y = 0.
// Triangulos lineares (a ordem de aproximacao e escolhida no NeoPZ: u P2, p_w e p_a P1).
// Gerar a malha: gmsh -2 coluna.geo -o coluna.msh
SetFactory("Built-in");

W = 0.10; // largura (m)
H = 1.00; // altura (m)
nx = 4;   // quadrados na largura
ny = 40;  // quadrados na altura

Point(1) = {0, 0, 0};
Point(2) = {W, 0, 0};
Point(3) = {W, H, 0};
Point(4) = {0, H, 0};

Line(1) = {1, 2}; // base
Line(2) = {2, 3}; // lateral direita
Line(3) = {3, 4}; // topo
Line(4) = {4, 1}; // lateral esquerda

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// malha estruturada: nx x ny quadrados, cada um dividido em dois triangulos pela diagonal
Transfinite Curve{1, 3} = nx + 1;
Transfinite Curve{2, 4} = ny + 1;
Transfinite Surface{1} = {1, 2, 3, 4} Left;

Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 4.1;

Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2};
Physical Curve("Topo", 4) = {3};
Physical Curve("Esquerda", 5) = {4};
