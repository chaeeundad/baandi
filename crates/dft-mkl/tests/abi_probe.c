/* Compile against user-installed headers; no MKL declarations are copied. */
#include <mkl.h>
#include <mkl_dfti.h>
#include <stdio.h>
_Static_assert(sizeof(MKL_LONG) == sizeof(long), "DFTI long ABI");
_Static_assert(sizeof(MKL_INT) == sizeof(int), "BLAS LP64 ABI");
_Static_assert(sizeof(MKL_Complex16) == 2*sizeof(double), "complex ABI");
_Static_assert(DFTI_DOUBLE == 36 && DFTI_COMPLEX == 32, "domain constants");
_Static_assert(DFTI_FORWARD_SCALE == 4 && DFTI_PLACEMENT == 11, "parameter constants");
_Static_assert(DFTI_INPLACE == 43 && DFTI_NO_ERROR == 0, "value constants");
_Static_assert(CblasColMajor == 102 && CblasNoTrans == 111 && CblasConjTrans == 113, "CBLAS enums");
int main(void) {
    printf("PASS MKL LP64 ABI constants; MKL_LONG=%zu MKL_INT=%zu complex=%zu bytes\n",
        sizeof(MKL_LONG), sizeof(MKL_INT), sizeof(MKL_Complex16));
    return 0;
}
