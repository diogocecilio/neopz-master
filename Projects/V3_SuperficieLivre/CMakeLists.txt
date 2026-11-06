# V3_SuperficieLivre: percolacao permanente com superficie livre (variante iii): barragem retangular com agua a jusante (saida 8,9 m acima do nivel de jusante) e barragem com face inclinada (face de percolacao de 2,4 m).
# Referencia: Sharma, Fujisawa & Murakami (2021), Space-Time Finite Element Method for Transient and Unconfined Seepage Flow Analysis, Finite Elem. Anal. Des. 197:103632.
add_verif_exe(V3_SuperficieLivre main.cpp problema.h)
add_verif_test(V3_testes testes.cpp problema.h)
