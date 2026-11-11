// Talude de de Souza Neto, Peric & Owen (Fig. 8.30): H = 10 m a 45 graus, dominio 75 m (40 m a esquerda, 30 m a
// direita), crista de 35 m. E a altura do exemplo de Cho da Fig. 19 de Vargas Ceron et al. (2025): Gamma = 1,777.
// Malha grossa como a ceron_h5_a45_g15 com o dobro do tamanho (H dobrou): hSup = 3 m perto da superficie.
SetFactory("OpenCASCADE");
hSup = 3.0; hMax = 8.0;
Point(1) = {0, 0, 0};    // base, esquerda
Point(2) = {75, 0, 0};   // base, direita
Point(3) = {75, 30, 0};  // terreno do pe, direita
Point(4) = {45, 30, 0};  // pe
Point(5) = {35, 40, 0};  // crista
Point(6) = {0, 40, 0};   // crista, esquerda
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
Field[2].DistMin = 2.0; Field[2].DistMax = 15.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6; Mesh.ElementOrder = 1; Mesh.MshFileVersion = 4.1;
Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2};
Physical Curve("Superficie", 4) = {3, 4, 5};
Physical Curve("Esquerda", 5) = {6};
