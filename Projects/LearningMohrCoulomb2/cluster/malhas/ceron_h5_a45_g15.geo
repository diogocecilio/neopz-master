// Malha inicial ceron_h5_a45_g15 (126 triangulos com o gmsh 4.15.2) do estudo de niveis de refinamento em volta da banda de
// ruptura (../../../SlopeRandomFieldRefina/niveis.py): a ceron_h5_a45_foco375 sem o foco (sem os pontos das bandas e os
// seus campos de tamanho), so com hSup = 1,5 m perto da superficie e hMax = 3 m. Mesmo contorno, mesmos 6 cantos e
// mesmos grupos fisicos da foco375: o retangulo do Karhunen-Loeve (0 a 25 m por 0 a 10 m) e o campo aleatorio de cada
// realizacao sao os mesmos nas duas malhas e nas malhas refinadas a partir delas (refina_banda.py).
// Talude de Vargas Ceron et al. (2025), Secao 6 (Fig. 24): H = 5 m, beta = 45 graus. Dominio como o dos exemplos da
// Secao 5.3 do artigo (Figs. 15 e 19): 10 m de crista antes da face e 10 m de terreno depois do pe, base 5 m abaixo do
// pe. Base em y = 0, pe em (15, 5), crista em (10, 10). Refinamento como o de rahardjo_a45.geo, na escala de H (a
// metade): hSup ate 1 m da superficie, crescendo ate hMax a 6 m. Grupos fisicos com os nomes usados em ReadGmshMesh.
SetFactory("OpenCASCADE");
hSup = 1.5; hMax = 3.0;
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
Field[1] = Distance; Field[1].CurvesList = {3, 4, 5}; Field[1].Sampling = 400;
Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = hSup; Field[2].SizeMax = hMax;
Field[2].DistMin = 1.0; Field[2].DistMax = 6.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6; Mesh.ElementOrder = 1; Mesh.MshFileVersion = 4.1;
Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2};
Physical Curve("Superficie", 4) = {3, 4, 5};
Physical Curve("Esquerda", 5) = {6};
