// Talude (de Souza Neto, Peric & Owen, Fig. 8.30-8.31): altura h = 10 m, inclinacao 45 graus,
// dominio 75 m x 40 m (esquerda) / 30 m (direita), crista de 35 m. Unidades: m. y para cima.
// Triangulos lineares (a ordem de aproximacao e escolhida no NeoPZ), refinados junto a face,
// a crista e ao pe, onde se forma a superficie de ruptura.
SetFactory("Built-in");

hFine = 0.8;   // tamanho junto a face do talude
hCoarse = 7.0; // tamanho longe do talude

Point(1) = {0, 0, 0, hCoarse};
Point(2) = {75, 0, 0, hCoarse};
Point(3) = {75, 30, 0, hCoarse};
Point(4) = {45, 30, 0, hFine};   // pe (C)
Point(5) = {35, 40, 0, hFine};   // crista (A)
Point(6) = {0, 40, 0, hCoarse};

Line(1) = {1, 2}; // base
Line(2) = {2, 3}; // lateral direita
Line(3) = {3, 4}; // terreno a jusante
Line(4) = {4, 5}; // face do talude
Line(5) = {5, 6}; // crista
Line(6) = {6, 1}; // lateral esquerda

Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// tamanho cresce com a distancia a face do talude (a superficie de ruptura passa por ela)
Field[1] = Distance;
Field[1].CurvesList = {4};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = hFine;
Field[2].SizeMax = hCoarse;
Field[2].DistMin = 1.0;
Field[2].DistMax = 25.0;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.Algorithm = 6;      // Frontal-Delaunay: triangulos de boa qualidade
Mesh.Smoothing = 10;
Mesh.ElementOrder = 1;
Mesh.MshFileVersion = 4.1;

Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2};
Physical Curve("Superficie", 4) = {3, 4, 5};
Physical Curve("Esquerda", 5) = {6};
