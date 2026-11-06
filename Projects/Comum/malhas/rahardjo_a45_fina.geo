// Talude homogeneo de Rahardjo et al. (2007), Fig. 2: Hs = 10 m, alpha = 45 graus, plato de 3Hs na crista e no pe,
// base em y = 0, pe em y = 2Hs, crista em y = 3Hs. Mesmo refinamento da malha "media" de talude3f/malha_talude.py:
// hSup ate 2 m de profundidade, crescendo ate hMax a 2 m + Hs. Grupos fisicos com os nomes usados em ReadGmshMesh.
SetFactory("OpenCASCADE");
hSup = 0.375; hMax = 1.5;
Point(1) = {0, 0, 0};    // g
Point(2) = {70, 0, 0};   // f
Point(3) = {70, 20, 0};  // d
Point(4) = {40, 20, 0};  // c (pe)
Point(5) = {30, 30, 0};  // b (crista)
Point(6) = {0, 30, 0};   // a
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
Field[2].DistMin = 2.0; Field[2].DistMax = 12.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6; Mesh.ElementOrder = 1; Mesh.MshFileVersion = 4.1;
Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2};
Physical Curve("Superficie", 4) = {3, 4, 5};
Physical Curve("Esquerda", 5) = {6};
