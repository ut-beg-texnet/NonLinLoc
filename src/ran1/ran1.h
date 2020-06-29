
#ifdef EXTERN_MODE
#define EXTERN_TXT extern
#else
#define EXTERN_TXT
#endif

EXTERN_TXT int RanSeed;



int get_rand_int(const int, const int);
double get_rand_double(const double, const double);
void test_rand_int();


/*/////////////////////////////////////////////////////////////////////////// */
/* Random number generator */

/*** function to get random double between xmin and xmax */

/* UNIX */
/*
#define RAND_MAX1 1.0
#define SRAND_FUNC(x) srand48((long) x)
#define RAND_FUNC() drand48()
 */

/* VC++ rand() */
/*#define SRAND_FUNC(x) srand((unsigned int) x)
#define RAND_FUNC() rand()
#define RAND_MAX1 RAND_MAX*/

/* Num Rec ran1() */
/*
#define SRAND_FUNC(x) seed_ran1((int) x)
#define RAND_FUNC() get_ran1()
#define RAND_MAX1 1.0
 */

/* UNI */

#define SRAND_FUNC(x) rinit((int) x)
#define RAND_FUNC() uni()
#define RAND_MAX1 1.0


double get_rand_double(const double xmin, const double xmax);
int get_rand_int(const int imin, const int imax);
void test_rand_int();

/* */
/*/////////////////////////////////////////////////////////////////////////// */






/*//////// UNI stuff */

double uni(void);
void rstart(int i, int j, int k, int l);
void rinit(int ijkl);
