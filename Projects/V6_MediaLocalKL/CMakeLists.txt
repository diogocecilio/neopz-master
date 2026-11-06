# V6_MediaLocalKL: o gerador de Karhunen-Loeve: o desvio da media do campo ao longo de um arco circular contra a funcao de variancia (Eqs. 3 a 15), independente dos elementos finitos.
# Referencia: Yuan, Hu, Mei, Wang & Wang (2021), Slope reliability analysis based on curvilinear local averaging of a 2-D random field, Comput. Geotech. 137:104247.
add_verif_exe(V6_MediaLocalKL main.cpp problema.h)
add_verif_test(V6_testes testes.cpp problema.h)
