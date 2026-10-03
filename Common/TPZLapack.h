#ifndef TPZLapackH
#define TPZLapackH

#ifdef USING_LAPACK
#ifdef MKLLAPACK
#include <mkl.h>
typedef MKL_Complex16 vardoublecomplex;
typedef MKL_Complex8 varfloatcomplex;

#elif MACOSX
#include <Accelerate/Accelerate.h>
typedef __CLPK_doublecomplex vardoublecomplex;
typedef __CLPK_complex varfloatcomplex;

#else
#include <complex>

#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>

typedef lapack_complex_double vardoublecomplex;
typedef lapack_complex_float varfloatcomplex;

#include "lapacke.h"
#include <cblas.h>
#endif

// LAPACK >= 3.9.1 (lapack.h com LAPACK_FORTRAN_STRLEN_END) declara os comprimentos ocultos das strings
// Fortran nas rotinas com argumentos character; as macros LAPACK_xxx do lapack.h os acrescentam.
// Use PZ_LAPACK(dgesvd)(...) no lugar de dgesvd_(...).
#if defined(LAPACK_FORTRAN_STRLEN_END) && !defined(MKLLAPACK)
#define PZ_LAPACK(name) LAPACK_##name
#else
#define PZ_LAPACK(name) name##_
#endif
#endif

#endif
