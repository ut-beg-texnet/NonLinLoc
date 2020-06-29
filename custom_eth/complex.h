#ifndef _NR_COMPLEX_H_
#define _NR_COMPLEX_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {float r,i;} fcomplex;
typedef struct DCOMPLEX {double r,i;} dcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

dcomplex Cadd();
dcomplex Csub();
dcomplex Cmul();
dcomplex Complex();
dcomplex Conjg();
dcomplex Cdiv();
double Cabs();
dcomplex Csqrt();
dcomplex RCmul();

#endif /* _NR_COMPLEX_H_ */
