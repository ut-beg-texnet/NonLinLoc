

#include <stdlib.h>
#include <stdio.h>
#include "ran1.h"




/*//////// Anthony Lomax stuff */


/*** function to get random integer between imin and imax */

 int get_rand_int(const int imin, const int imax)
{

	return( imin + (int) ( (double) RAND_FUNC() * (double) (imax - imin + 1)
			/ (double) RAND_MAX1 ) );
}


/*** function to get random double between xmin and xmax */

 double get_rand_double(const double xmin, const double xmax)
{

	return( xmin + (double) RAND_FUNC()
		* (xmax - xmin) / (double) RAND_MAX1 );
}





/*** function to test random number generator */

#define NUM_BIN 16

void test_rand_int()
{
	long imax = 15, nmax = 32000;
	long n, m, itest;
	long ibin[NUM_BIN], ibinmax[NUM_BIN];

	for (n = 0; n < NUM_BIN; n++) {
		ibin[n] = 0;
		ibinmax[n] = (n + 1) * ((imax + 1) / NUM_BIN) - 1;
	}


	for (n = 0; n < nmax; n++) {
		itest = get_rand_int(0, imax);
		for (m = -1; itest > ibinmax[++m];);
		ibin[m]++;
	}

	fprintf(stdout,
		"\nRandom function test (val= 0 - %ld, samples= %ld)\n", imax, nmax);
	fprintf(stdout,
		"  RAND_MAX1 is %ld (%le)\n",
			(long) RAND_MAX1, (double) RAND_MAX1);
	fprintf(stdout, "  Bin 0-%ld  N=%ld\n", ibinmax[0], ibin[0]);
	for (n = 1; n < NUM_BIN; n++) {
		fprintf(stdout,
			"  Bin %ld-%ld  N=%ld\n", ibinmax[n - 1] + 1, ibinmax[n],
			ibin[n]);
	}

}








/*//////// UNI stuff */



/*
 *	C version of Marsaglia's UNI random number generator
 *	More or less transliterated from the Fortran -- with 1 bug fix
 *	Hence horrible style
 *
 *	Features:
 *		ANSI C
 *		not callable from Fortran (yet)
 */

/*	As far as I can tell, use init(ijkl) to initialise,
 *	then uni() to generate a random number in the range (0-1).
 */

/*
 *	Global variables for rstart & uni
 */

double uni_u[98];	/* Was U(97) in Fortran version -- too lazy to fix */
double uni_c, uni_cd, uni_cm;
int uni_ui, uni_uj;

 double uni(void)
{
	double luni;			/* local variable for uni */

	luni = uni_u[uni_ui] - uni_u[uni_uj];
	if (luni < 0.0)
		luni += 1.0;
	uni_u[uni_ui] = luni;
	if (--uni_ui == 0)
		uni_ui = 97;
	if (--uni_uj == 0)
		uni_uj = 97;
	if ((uni_c -= uni_cd) < 0.0)
		uni_c += uni_cm;
	if ((luni -= uni_c) < 0.0)
		luni += 1.0;
	return (double) luni;
}

 void rstart(int i, int j, int k, int l)
{
	int ii, jj, m;
	double s, t;

	for (ii = 1; ii <= 97; ii++) {
		s = 0.0;
		t = 0.5;
		for (jj = 1; jj <= 24; jj++) {
			m = ((i*j % 179) * k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53*l+1) % 169;
			if (l*m % 64 >= 32)
				s += t;
			t *= 0.5;
		}
		uni_u[ii] = s;
	}
	uni_c  = 362436.0   / 16777216.0;
	uni_cd = 7654321.0  / 16777216.0;
	uni_cm = 16777213.0 / 16777216.0;
	uni_ui = 97;	/*  There is a bug in the original Fortran version */
	uni_uj = 33;	/*  of UNI -- i and j should be SAVEd in UNI()     */
}


/* ~rinit: this takes a single integer in the range
		0 <= ijkl <= 900 000 000
	and produces the four smaller integers needed for rstart. It is
 *	based on the ideas contained in the RMARIN subroutine in
 *		F. James, "A Review of Pseudorandom Number Generators",
 *			Comp. Phys. Commun. Oct 1990, p.340
 *	To reduce the modifications to the existing code, rinit now
 *	takes the role of a preprocessor for rstart.
 *
 *	This is useful for the parallel version of the code as James
 *	states that any integer ijkl will produce a statistically
 *	independent sequence of random numbers.
 *
 *     Very funny. If that statement was worth anything he would have provided
 *     a proof to go with it. spb 12/12/90
 */

void rinit(int ijkl)
{
	int i, j, k, l, ij, kl;

	/* check ijkl is within range */
	if( (ijkl < 0) || (ijkl > 900000000) )
		{
		printf("rinit: ijkl = %d -- out of range\n\n", ijkl);
		exit(3);
               	}

/*        printf("rinit: seed_ijkl = %d\n", ijkl); */

	/* decompose the long integer into the the equivalent four
	 * integers for rstart. This should be a 1-1 mapping
 	 *	ijkl <--> (i, j, k, l)
	 * though not quite all of the possible sets of (i, j, k, l)
	 * can be produced.
	 */

	ij = ijkl/30082;
	kl = ijkl - (30082 * ij);

	i = ((ij/177) % 177) + 2;
	j = (ij % 177) + 2;
	k = ((kl/169) % 178) + 1;
	l = kl % 169;

	if( (i <= 0) || (i > 178) )
		{
		printf("rinit: i = %d -- out of range\n\n", i);
		exit(3);
               	}

	if( (j <= 0) || (j > 178) )
		{
		printf("rinit: j = %d -- out of range\n\n", j);
		exit(3);
               	}

	if( (k <= 0) || (k > 178) )
		{
		printf("rinit: k = %d -- out of range\n\n", k);
		exit(3);
               	}

	if( (l < 0) || (l > 168) )
		{
		printf("rinit: l = %d -- out of range\n\n", l);
		exit(3);
               	}

	if (i == 1 && j == 1 && k == 1)
		{
                printf("rinit: 1 1 1 not allowed for 1st 3 seeds\n\n");
		exit(4);
                }

/*        printf("rinit: initialising RNG via rstart(%d, %d, %d, %d)\n",
				i, j, k, l); */

        rstart(i, j, k, l);

}


