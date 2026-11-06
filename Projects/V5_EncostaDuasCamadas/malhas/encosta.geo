// V5_EncostaDuasCamadas: encosta em duas camadas de Abbasov et al. (2024), Fig. 2 (benchmark de Lu et al. 2012 /
// Shao et al. 2015). Dominio 42 m x 25 m; crista z = 25 m ate x = 19 m, talude de 6 m em 15 m (21,8 graus) ate o pe
// em x = 34 m, z = 19 m; camada superior de 2 m (espessura vertical) acima da interface (a superficie 2 m abaixo).
// Zona de investigacao: 15 <= x <= 37 m, z >= 17 m (22 m de largura, 8 m sob a crista, 2 m sob o pe; 113 m2).
// As arestas da interface e da zona de investigacao sao arestas da malha (as integrais na zona sao exatas).
// Camada superior: triangulos estruturados (transfinita), nz linhas na vertical; camada inferior: nao estruturada.
// Gerar: gmsh -2 encosta.geo -format msh41 -o encosta.msh [-setnumber nz 16 -setnumber hx 0.25 ...]
DefineConstant[ nz = {16, Name "linhas na camada superior"} ];
DefineConstant[ hx = {0.25, Name "espacamento horizontal na zona (m)"} ];
DefineConstant[ hxe = {0.75, Name "espacamento horizontal a esquerda da zona (m)"} ];
DefineConstant[ hxd = {0.5, Name "espacamento horizontal a direita da zona (m)"} ];
DefineConstant[ hz = {0.3, Name "tamanho na camada inferior dentro da zona (m)"} ];
DefineConstant[ hb = {3.0, Name "tamanho na base (m)"} ];

// pontos (x, z)
Point(1) = {0, 0, 0, hb};
Point(2) = {42, 0, 0, hb};
Point(3) = {42, 17, 0, hxd};   // lateral direita: carga hidrostatica abaixo, face de percolacao acima
Point(4) = {42, 19, 0, hxd};
Point(5) = {37, 19, 0, hx};    // borda direita da zona na superficie
Point(6) = {34, 19, 0, hx};    // pe do talude
Point(7) = {19, 25, 0, hx};    // borda da crista
Point(8) = {15, 25, 0, hx};    // borda esquerda da zona na superficie
Point(9) = {0, 25, 0, hxe};
Point(10) = {0, 23, 0, hxe};   // interface na lateral esquerda
Point(11) = {15, 23, 0, hx};
Point(12) = {19, 23, 0, hx};
Point(13) = {34, 17, 0, hx};
Point(14) = {37, 17, 0, hx};
Point(15) = {15, 17, 0, hz};   // canto inferior esquerdo da zona

// contornos externos
Line(1) = {1, 2};              // base
Line(2) = {2, 3};              // direita, z <= 17
Line(3) = {3, 4};              // direita, 17 <= z <= 19
Line(4) = {4, 5};              // superficie (pe, fora da zona)
Line(5) = {5, 6};              // superficie (pe, na zona)
Line(6) = {6, 7};              // superficie (talude)
Line(7) = {7, 8};              // superficie (crista, na zona)
Line(8) = {8, 9};              // superficie (crista, fora da zona)
Line(9) = {9, 10};             // esquerda, camada superior
Line(10) = {10, 1};            // esquerda, camada inferior
// interface
Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 3};
// verticais da camada superior (limites dos blocos transfinitos e da zona)
Line(16) = {11, 8};
Line(17) = {12, 7};
Line(18) = {13, 6};
Line(19) = {14, 5};
// limites da zona na camada inferior
Line(20) = {15, 11};
Line(21) = {15, 13};

// camada superior: cinco quadrilateros
Curve Loop(1) = {11, 16, 8, 9};        Plane Surface(1) = {1};
Curve Loop(2) = {12, 17, 7, -16};      Plane Surface(2) = {2};
Curve Loop(3) = {13, 18, 6, -17};      Plane Surface(3) = {3};
Curve Loop(4) = {14, 19, 5, -18};      Plane Surface(4) = {4};
Curve Loop(5) = {15, 3, 4, -19};       Plane Surface(5) = {5};
// camada inferior: a parte na zona e o resto
Curve Loop(6) = {21, -13, -12, -20};   Plane Surface(6) = {6};
Curve Loop(7) = {1, 2, -15, -14, -21, 20, -11, 10}; Plane Surface(7) = {7};

// divisoes (numero de nos = divisoes + 1)
Transfinite Curve{9, 16, 17, 18, 19, 3} = nz + 1;
Transfinite Curve{11, 8} = Ceil(15 / hxe) + 1;
Transfinite Curve{12, 7} = Ceil(4 / hx) + 1;
Transfinite Curve{13, 6} = Ceil(Sqrt(15^2 + 6^2) / hx) + 1;
Transfinite Curve{14, 5} = Ceil(3 / hx) + 1;
Transfinite Curve{15, 4} = Ceil(5 / hxd) + 1;
Transfinite Surface{1} = {10, 11, 8, 9} Alternate;
Transfinite Surface{2} = {11, 12, 7, 8} Alternate;
Transfinite Surface{3} = {12, 13, 6, 7} Alternate;
Transfinite Surface{4} = {13, 14, 5, 6} Alternate;
Transfinite Surface{5} = {14, 3, 4, 5} Alternate;

Physical Surface("Superior") = {1, 2, 3, 4, 5};
Physical Surface("Inferior") = {6, 7};
Physical Curve("Base") = {1};
Physical Curve("Direita") = {2, 3};
Physical Curve("Superficie") = {4, 5, 6, 7, 8};
Physical Curve("Esquerda") = {9, 10};
Mesh.Algorithm = 6;
