// Talude de Liu et al. (2024), Fig. 1: altura 5 m, inclinacao 1:2, comprimento 30 m (o mesmo de Griffiths e Fenton 2004,
// Cho 2010 e Jiang et al. 2014). Pontos (m): C(0,0), D(30,0), E(30,1) (nivel d'agua a direita), F(30,5), G(20,5) (pe),
// H(10,10) (crista), A(0,10), B(0,6) (nivel d'agua a esquerda). Triangulos de tamanho h (o artigo usa quadrilateros de
// 0,5 m x 0,5 m no SEEP/W). Uso: gmsh -2 talude_liu.geo -setnumber h 0.5 -format msh41 -o talude_liu_h05.msh
If (!Exists(h))
  h = 0.5;
EndIf
Point(1) = {0, 0, 0, h};   // C
Point(2) = {30, 0, 0, h};  // D
Point(3) = {30, 1, 0, h};  // E
Point(4) = {30, 5, 0, h};  // F
Point(5) = {20, 5, 0, h};  // G (pe)
Point(6) = {10, 10, 0, h}; // H (crista)
Point(7) = {0, 10, 0, h};  // A
Point(8) = {0, 6, 0, h};   // B
Line(1) = {1, 2}; // base CD
Line(2) = {2, 3}; // direita DE (carga 1 m)
Line(3) = {3, 4}; // direita EF (impermeavel)
Line(4) = {4, 5}; // terreno do pe FG
Line(5) = {5, 6}; // face GH
Line(6) = {6, 7}; // crista HA
Line(7) = {7, 8}; // esquerda AB (impermeavel)
Line(8) = {8, 1}; // esquerda BC (carga 6 m)
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
Mesh.Algorithm = 6; Mesh.ElementOrder = 1; Mesh.MshFileVersion = 4.1;
Mesh.MeshSizeExtendFromBoundary = 1; Mesh.MeshSizeFromCurvature = 0;
Physical Surface("Solo", 1) = {1};
Physical Curve("Base", 2) = {1};
Physical Curve("Direita", 3) = {2, 3};
Physical Curve("Superficie", 4) = {4, 5, 6};
Physical Curve("Esquerda", 5) = {7, 8};
