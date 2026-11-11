// Talude do relatorio (H = 5 m, beta = 45 graus) num dominio de 10H: 50 m de crista atras da face, 50 m de terreno
// depois do pe e base 50 m abaixo do pe. Cobre o raio R_e = sqrt(H^2 + (L_m + H/tan(beta))^2) ~ 55 m do campo de
// percolacao do artigo (L_m = 10H). Perto da face como a ceron_h5_a45_g15 (1,5 m), crescendo ate 10 m longe dela.
SetFactory("OpenCASCADE");
hSup = 1.5; hMax = 10.0;
Point(1) = {0, 0, 0};     // base, esquerda
Point(2) = {105, 0, 0};   // base, direita
Point(3) = {105, 50, 0};  // terreno do pe, direita
Point(4) = {55, 50, 0};   // pe
Point(5) = {50, 55, 0};   // crista
Point(6) = {0, 55, 0};    // crista, esquerda
Line(1) = {1, 2};  // base
Line(2) = {2, 3};  // lateral direita
Line(3) = {3, 4};  // terreno do pe
Line(4) = {4, 5};  // face
Line(5) = {5, 6};  // crista
Line(6) = {6, 1};  // lateral esquerda
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
Field[1] = Distance; Field[1].CurvesList = {4}; Field[1].Sampling = 400;
Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = hSup; Field[2].SizeMax = hMax;
Field[2].DistMin = 2.0; Field[2].DistMax = 40.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6; Mesh.ElementOrder = 1; Mesh.MshFileVersion = 4.1;
Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2};
Physical Curve("Superficie", 4) = {3, 4, 5};
Physical Curve("Esquerda", 5) = {6};
