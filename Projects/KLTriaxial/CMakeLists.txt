# Sect. 7.1 (v1.0): triaxial tests of Krabbenhoft and Lyamin (2012) with one Hex20-Hex8 element, drained and
# undrained, NC and OCR = 10, against the closed forms of Appendix B
add_mcc_example(KLTriaxial main.cpp KLTriaxial.h)
# IntegrationSchemes.h (header-only) provides the clay, the undrained closed form and Brent's method
target_include_directories(KLTriaxial PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../IntegrationSchemes)
