// Talude de Vargas Ceron et al. (2025), Secao 6 (H = 5 m, beta = 45 graus): o mesmo poligono e os mesmos grupos fisicos
// de ../../LearningMohrCoulomb2/cluster/malhas/ceron_h5_a45.geo, com triangulos de tamanho uniforme h = 0.5 m em todo o
// dominio (sem refinamento junto a superficie). Base em y = 0, pe em (15, 5), crista em (10, 10); area 187,5 m2.
// Gerada por gera_malhas.py (modulo gmsh do Python, 4.15, uma thread).
SetFactory("OpenCASCADE");
h = 0.5;
Point(1) = {0, 0, 0};    // base, esquerda
Point(2) = {25, 0, 0};   // base, direita
Point(3) = {25, 5, 0};   // terreno do pe, direita
Point(4) = {15, 5, 0};   // pe
Point(5) = {10, 10, 0};  // crista
Point(6) = {0, 10, 0};   // crista, esquerda
Line(1) = {1, 2};  // base
Line(2) = {2, 3};  // lateral direita
Line(3) = {3, 4};  // terreno do pe
Line(4) = {4, 5};  // face
Line(5) = {5, 6};  // crista
Line(6) = {6, 1};  // lateral esquerda
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
Mesh.MeshSizeMin = h; Mesh.MeshSizeMax = h;
Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0; Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.Algorithm = 6; Mesh.ElementOrder = 1; Mesh.MshFileVersion = 4.1;
Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2};
Physical Curve("Superficie", 4) = {3, 4, 5};
Physical Curve("Esquerda", 5) = {6};
