# V5_EncostaDuasCamadas: encosta em duas camadas (benchmark de Lu et al. 2012): 10 anos a 600 mm/ano e 20 mm/h por 15 h, chi = S_e, mapas de umidade e do fator de seguranca local (LFS).
# Referencia: Abbasov et al. (2024), Modeling rainfall-induced landslide using the concept of local factor of safety: uncertainty propagation and sensitivity analysis, Comput. Geotech. 167:106102.
add_verif_exe(V5_EncostaDuasCamadas main.cpp problema.h)
add_verif_test(V5_testes testes.cpp problema.h)
