// Talude do Ceron (H = 5 m, 45 graus, crista em (10, 10), pe em (15, 5)) num dominio grande, para a zona de fluxo do
// artigo (R_e com L_m = 10 H, cerca de 55 m em torno da crista): 50 m de crista antes da face, 50 m de terreno depois do
// pe, base 50 m abaixo da crista. hFoco = 0,25 m na regiao das bandas de ruptura (caixa x 7,4 a 17,5, y 3,5 a 10,1),
// hSup = 0,6 m perto da superficie, crescendo ate hMax = 6 m a 25 m dela.
SetFactory("OpenCASCADE");
hFoco = 0.25; hSup = 0.6; hMax = 6.0;
Point(1) = {-40, -45, 0};
Point(2) = {65, -45, 0};
Point(3) = {65, 5, 0};
Point(4) = {15, 5, 0};
Point(5) = {10, 10, 0};
Point(6) = {-40, 10, 0};
Line(1) = {1, 2};  // base
Line(2) = {2, 3};  // lateral direita
Line(3) = {3, 4};  // terreno do pe
Line(4) = {4, 5};  // face
Line(5) = {5, 6};  // crista
Line(6) = {6, 1};  // lateral esquerda
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
Field[1] = Box; Field[1].VIn = hFoco; Field[1].VOut = hMax;
Field[1].XMin = 7.4; Field[1].XMax = 17.5; Field[1].YMin = 3.5; Field[1].YMax = 10.1; Field[1].Thickness = 3;
Field[2] = Distance; Field[2].CurvesList = {3, 4, 5}; Field[2].Sampling = 4000;
Field[3] = Threshold; Field[3].InField = 2; Field[3].SizeMin = hSup; Field[3].SizeMax = hMax;
Field[3].DistMin = 1.0; Field[3].DistMax = 25.0;
Field[4] = Min; Field[4].FieldsList = {1, 3};
Background Field = 4;
Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6; Mesh.ElementOrder = 1; Mesh.MshFileVersion = 4.1;
Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2};
Physical Curve("Superficie", 4) = {3, 4, 5};
Physical Curve("Esquerda", 5) = {6};
