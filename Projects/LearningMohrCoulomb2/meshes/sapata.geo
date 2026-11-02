// Sapata corrida (Fig. 4 do artigo de Lira Cecilio): dominio 5 m x 5 m, meia sapata de 0,5 m
// no canto superior esquerdo (simetria em x = 0). Unidades: m. y para cima.
// Triangulos lineares (a ordem de aproximacao e escolhida no NeoPZ), refinados sob a sapata
// e na borda dela (singularidade em x = 0,5 m).
SetFactory("Built-in");

hFine = 0.04;  // tamanho na borda da sapata
hCoarse = 0.8; // tamanho longe da sapata

Point(1) = {0, 0, 0, hCoarse};
Point(2) = {5, 0, 0, hCoarse};
Point(3) = {5, 5, 0, hCoarse};
Point(4) = {0.5, 5, 0, hFine}; // borda da sapata
Point(5) = {0, 5, 0, hFine};   // eixo de simetria

Line(1) = {1, 2}; // base
Line(2) = {2, 3}; // lateral direita
Line(3) = {3, 4}; // superficie livre
Line(4) = {4, 5}; // sapata
Line(5) = {5, 1}; // eixo de simetria

Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

// tamanho cresce com a distancia a sapata e, mais devagar, com a distancia a sua borda
Field[1] = Distance;
Field[1].CurvesList = {4};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 2 * hFine;
Field[2].SizeMax = hCoarse;
Field[2].DistMin = 0.1;
Field[2].DistMax = 3.0;
Field[3] = Distance;
Field[3].PointsList = {4};
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = hFine;
Field[4].SizeMax = hCoarse;
Field[4].DistMin = 0.05;
Field[4].DistMax = 2.5;
Field[5] = Min;
Field[5].FieldsList = {2, 4};
Background Field = 5;
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
Physical Curve("Superficie", 4) = {3};
Physical Curve("Simetria", 5) = {5};
Physical Curve("Sapata", 6) = {4};
