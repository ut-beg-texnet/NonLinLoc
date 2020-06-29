
// Linux fix
#define xmalloc malloc
#define xrealloc realloc


#include <stdio.h>
#include <math.h>

#define bool int
#define false 0
#define true 1
#define progname "do_spline"
#define FUZZ 0.0000001		/* potential roundoff error */
#define TRIG_ARG_MIN 0.001
#define EXIT_FAILURE -1

inline bool is_monotonic (int n, double *t);
inline double interpolate (int n, double *t, double *y, double *z, double x, double tension, bool periodic);
inline double quotient_sinh_func (double x, double y);
inline double quotient_sin_func (double x, double y);
inline void non_monotonic_error (void);
inline void fit (int n, double *t, double *y, double *z, double k, double tension, bool periodic);
inline double sin_func (double x);
inline double sinh_func (double x);
inline double tan_func (double x);
inline double tanh_func (double x);



/* do_spline() is the main routine for piecewise cubic spline
   interpolation, supporting both periodicity and a user-specified boundary
   condition parameter.  Nonzero tension may be specified, in which case
   the interpolate() routine, which this calls, will use not cubic
   polynomials but rather expressions involving hyperbolic sines.

   t[] and y[] are the arrays in which the abscissa and ordinate values of
   the user-specified data points are stored, and z[] is the array in which
   the 2nd derivatives at the knots (data points in the interior of the
   interval) will be stored.  used+1 is the effective size of each of these
   arrays.  The number of points supplied by the user was used+1 in the
   non-periodic case.  It was used+0 in the periodic case.  

   The reason that the number of elements is greater by one in the periodic
   case is that the first user-supplied data point occurs also at the end.
   In fact, in the periodic case this function will increment the size of
   the array once more, since the periodic interpolation algorithm requires
   the first two data points, not just the first, to appear at the end. */

inline void
do_spline (int used, int len, double **t, int ydimension, double **y, double **z, 
	   double tension, bool periodic, bool spec_boundary_condition,
	   double k, int precision, double first_t, double last_t, 
	   double spacing_t, int no_of_intervals, bool spec_first_t, 
	   bool spec_last_t, bool spec_spacing_t, 
	   bool spec_no_of_intervals, bool suppress_abscissa,
/* !!!AJL */
		double *do_spline_x, double **do_spline_yy
/* !!!AJL */
		)
{
  int range_count = 0;		/* number of req'd datapoints out of range */
  int lastval = 0;		/* last req'd point = 1st/last data point? */
  int i;

/* !!!AJL */
/*
printf("INTO do_spline: %d used, %d len, double **t, %d ydimension, double **y, double **z, 
	   %lf tension, bool %d periodic, bool %d spec_boundary_condition,
	   double k, %d precision, double first_t, double last_t, 
	   %lf spacing_t, %d no_of_intervals, bool %d spec_first_t, 
	   bool %d spec_last_t, bool %d spec_spacing_t, 
	   bool %d spec_no_of_intervals, bool %d suppress_abscissa\n",
	used, len, ydimension, tension, periodic, spec_boundary_condition,
	   precision, spacing_t, no_of_intervals, spec_first_t, 
	   spec_last_t, spec_spacing_t, 
	   spec_no_of_intervals, suppress_abscissa);
*/
/* !!!AJL */


  if (used + 1 == 0)		/* zero data points in array */
    /* don't output anything (i.e. effectively output a null dataset) */
    return;

  if (used+1 == 1)		/* a single data point in array */
    {
      fprintf (stderr, 
	       "%s: cannot construct a spline from a single data point\n", 
	       progname);
      /* don't output anything (i.e. effectively output a null dataset) */
      return;
    }

  if (!periodic && used+1 <= 2)
    {
      if (spec_boundary_condition)
	fprintf (stderr, 
		 "%s: only 2 data points, so ignoring specified boundary condition\n", 
		 progname);
      k = 0.0;
    }

  if (!is_monotonic (used, *t))
    non_monotonic_error();	/* self-explanatory */

  if (periodic)
    {
      bool print_warning = false;
      
      for (i = 0; i < ydimension; i++)
	{
	  if (y[i][used] != y[i][0])
	    print_warning = true;
	  y[i][used] = y[i][0];
	}
      if (print_warning == true)
	fprintf (stderr, "%s: setting final y value equal to initial to ensure periodicity\n", 
		 progname); 

      /* add pseudo-point at end (to accord with periodicity) */
      if (used + 1 >= len)
	{
	  len++;
	  *t = (double *)xrealloc (*t, sizeof(double) * len);
	  for (i = 0; i < ydimension; i++)
	    {
	      y[i] = (double *)xrealloc (y[i], sizeof(double) * len);
	      z[i] = (double *)xrealloc (z[i], sizeof(double) * len);
	    }
	}
      (*t)[used + 1] = (*t)[used] + ((*t)[1] - (*t)[0]);
      for (i = 0; i < ydimension; i++)
	y[i][used + 1] = y[i][1];
    }

  /* compute z[], array of 2nd derivatives at each knot */
  for (i = 0; i < ydimension; i++)
    fit (used + (periodic ? 1 : 0), /* include pseudo-point if any */
	 *t, y[i], z[i], k, tension, periodic);

  if (!spec_first_t) 
    first_t = (*t)[0];
  if (!spec_last_t) 
    last_t = (*t)[used];	/* used+1 data points in all */
  if (!spec_spacing_t) 
    {
      if (no_of_intervals > 0)
	spacing_t = (last_t - first_t) / no_of_intervals;
      else
	spacing_t = 0;		/* won't happen */
    }
  else				/* user specified spacing */
    {
      if ((last_t - first_t) * spacing_t < 0.0)
	{
	  fprintf (stderr, "%s: specified spacing is of wrong sign, corrected\n",
		   progname);
	  spacing_t = -spacing_t;
	}
      if (spec_no_of_intervals)
	fprintf (stderr, "%s: ignoring specified number of intervals\n",
		 progname);
      no_of_intervals = (int)(fabs((last_t - first_t) / spacing_t) + FUZZ);
    }

  if (last_t == (*t)[0])
    lastval = 1;
  else if (last_t == (*t)[used])
    lastval = 2;

  for (i = 0; i <= no_of_intervals; ++i)
    {
      double x;

      x = first_t + spacing_t * i;

      if (i == no_of_intervals)
	{
	  /* avoid numerical fuzz */
	  if (lastval == 1)	/* left end of input */
	    x = (*t)[0];
	  else if (lastval == 2) /* right end of input */
	    x = (*t)[used];
	}

      if (periodic || (x - (*t)[0]) * (x - (*t)[used]) <= 0)
	{
/* !!!AJL */
	  int j;
	  //double *yy;

	  //yy = (double *)xmalloc (sizeof(double) * ydimension); 
	  for (j = 0; j < ydimension; j++) {
	    do_spline_yy[j][i] = interpolate (used, *t, y[j], z[j], x, 
				 tension, periodic);
	    do_spline_x[i] = x; 
		}
	  //write_point (x, yy, ydimension, precision, suppress_abscissa);
	  //free (yy);
/* !!!AJL */
	}
      else
	range_count++;
    }

  switch (range_count)
    {
    case 0:
      break;
    case 1:
      fprintf (stderr, 
	       "%s: one requested point could not be computed (out of data range)\n", 
	       progname);
      break;
    default:
      fprintf (stderr, 
	       "%s: %d requested points could not be computed (out of data range)\n", 
	       progname, range_count);
      break;
    }
}




/* fit() computes the array z[] of second derivatives at the knots, i.e.,
   internal data points.  The abscissa array t[] and the ordinate array y[]
   are specified.  On entry, have n+1 >= 2 points in the t, y, z arrays,
   numbered 0..n.  The knots are numbered 1..n-1 as in Kincaid and Cheney.
   In the periodic case, the final knot, i.e., (t[n-1],y[n-1]), has the
   property that y[n-1]=y[0]; moreover, y[n]=y[1].  The number of points
   supplied by the user was n+1 in the non-periodic case, and n in the
   periodic case.  When this function is called, n>=1 in the non-periodic
   case, and n>=2 in the periodic case. */

/* Algorithm: the n-1 by n-1 tridiagonal matrix equation for the vector of
   2nd derivatives at the knots is reduced to upper diagonal form.  At that
   point the diagonal entries (pivots) of the upper diagonal matrix are in
   the vector u[], and the vector on the right-hand side is v[].  That is,
   the equation is of the form Ay'' = v, where a_(ii) = u[i], and a_(i,i+1)
   = alpha[i].  Here i=1..n-1 indexes the set of knots.  The matrix
   equation is solved by back-substitution for y''[], i.e., for z[]. */

void
#ifdef _HAVE_PROTOS
fit (int n, double *t, double *y, double *z, double k, double tension,
     bool periodic)
#else
fit (n, t, y, z, k, tension, periodic)
     int n;
     double *t, *y, *z;
     double k;			/* y''_1 = k y''_0, etc. */
     double tension;
     bool periodic;
#endif
{
  double *h, *b, *u, *v, *alpha, *beta;
  double *uu = NULL, *vv = NULL, *s = NULL;
  int i;

  if (n == 1)			/* exactly 2 points, use straight line */
    {
      z[0] = z[1] = 0.0;
      return;
    }

  h = (double *)xmalloc (sizeof(double) * n);
  b = (double *)xmalloc (sizeof(double) * n);
  u = (double *)xmalloc (sizeof(double) * n);
  v = (double *)xmalloc (sizeof(double) * n);
  alpha = (double *)xmalloc (sizeof(double) * n);
  beta = (double *)xmalloc (sizeof(double) * n);

  if (periodic)
    {
      s = (double *)xmalloc (sizeof(double) * n); 
      uu = (double *)xmalloc (sizeof(double) * n); 
      vv = (double *)xmalloc (sizeof(double) * n); 
    }

  for (i = 0; i <= n - 1 ; ++i)
    {
      h[i] = t[i + 1] - t[i];
      b[i] = 6.0 * (y[i + 1] - y[i]) / h[i]; /* for computing RHS */
    }

  if (tension < 0.0)		/* must rule out sin(tension * h[i]) = 0 */
    {
      for (i = 0; i <= n - 1 ; ++i)
	if (sin (tension * h[i]) == 0.0)
	  {
	    fprintf (stderr, "%s: error: specified negative tension value is singular\n", progname);
	    exit (EXIT_FAILURE);
	  }
    }
  if (tension == 0.0)
    {
      for (i = 0; i <= n - 1 ; ++i)
	{
	  alpha[i] = h[i];	/* off-diagonal = alpha[i] to right */
	  beta[i] = 2.0 * h[i];	/* diagonal = beta[i-1] + beta[i] */
	}
    }
  else
    if (tension > 0.0)
      /* `positive' (really real) tension, use hyperbolic trig funcs */
      {
	for (i = 0; i <= n - 1 ; ++i)
	  {
	    if (fabs(tension * h[i]) < TRIG_ARG_MIN)
	      /* hand-compute (6/x^2)(1-x/sinh(x)) and (3/x^2)(x/tanh(x)-1)
                 to improve accuracy */
	      {
		alpha[i] = h[i] * sinh_func (tension * h[i]);
		beta[i] = 2.0 * h[i] * tanh_func (tension * h[i]);
	      }
	    else
	      {
		alpha[i] = ((6.0 / (tension * tension))
			   * ((1.0 / h[i]) - tension / sinh (tension * h[i])));
		beta[i] = ((6.0 / (tension * tension))
			   * (tension / tanh (tension * h[i]) - (1.0 / h[i])));
	      }
	  }
      }
    else				/* tension < 0 */
      /* `negative' (really imaginary) tension,  use circular trig funcs */
      {
	for (i = 0; i <= n - 1 ; ++i)
	  {
	    if (fabs(tension * h[i]) < TRIG_ARG_MIN)
	      /* hand-compute (6/x^2)(1-x/sin(x)) and (3/x^2)(x/tan(x)-1)
                 to improve accuracy */
	      {
		alpha[i] = h[i] * sin_func (tension * h[i]);
		beta[i] = 2.0 * h[i] * tan_func (tension * h[i]);
	      }
	    else
	      {
		alpha[i] = ((6.0 / (tension * tension))
		           * ((1.0 / h[i]) - tension / sin (tension * h[i])));
		beta[i] = ((6.0 / (tension * tension))
			   * (tension / tan (tension * h[i]) - (1.0 / h[i])));
	      }
	  }
      }
  
  if (!periodic && n == 2)
      u[1] = beta[0] + beta[1] + 2 * k * alpha[0];
  else
    u[1] = beta[0] + beta[1] + k * alpha[0];

  v[1] = b[1] - b[0];
  
  if (u[1] == 0.0)
    {
      fprintf (stderr, 
	       "%s: error: as posed, problem of computing spline is singular\n",
	       progname);
      exit (EXIT_FAILURE);
    }

  if (periodic)
    {
      s[1] = alpha[0];
      uu[1] = 0.0;
      vv[1] = 0.0;
    }

  for (i = 2; i <= n - 1 ; ++i)
    {
      u[i] = (beta[i] + beta[i - 1]
	      - alpha[i - 1] * alpha[i - 1] / u[i - 1]
	      + (i == n - 1 ? k * alpha[n - 1] : 0.0));

      if (u[i] == 0.0)
	{
	  fprintf (stderr, 
		   "%s: error: as posed, problem of computing spline is singular\n",
		   progname);
	  exit (EXIT_FAILURE);
	}


      v[i] = b[i] - b[i - 1] - alpha[i - 1] * v[i - 1] / u[i - 1];

      if (periodic)
	{
	  s[i] = - s[i-1] * alpha[i-1] / u[i-1];
	  uu[i] = uu[i-1] - s[i-1] * s[i-1] / u[i-1];
	  vv[i] = vv[i-1] - v[i-1] * s[i-1] / u[i-1];
	}
    }
      
  if (!periodic)
    {
      /* fill in 2nd derivative array */
      z[n] = 0.0;
      for (i = n - 1; i >= 1; --i)
	z[i] = (v[i] - alpha[i] * z[i + 1]) / u[i];
      z[0] = 0.0;
      
      /* modify to include boundary condition */
      z[0] = k * z[1];
      z[n] = k * z[n - 1];
    }
  else		/* periodic */
    {
      z[n-1] = (v[n-1] + vv[n-1]) / (u[n-1] + uu[n-1] + 2 * s[n-1]);
      for (i = n - 2; i >= 1; --i)
	z[i] = ((v[i] - alpha[i] * z[i + 1]) - s[i] * z[n-1]) / u[i];

      z[0] = z[n-1];
      z[n] = z[1];
    }

  if (periodic)
    {
      free (vv);
      free (uu);
      free (s);
    }
  free (beta);
  free (alpha);
  free (v);
  free (u);
  free (b);
  free (h);
}

 

 
/* interpolate() computes an approximate ordinate value for a given
   abscissa value, given an array of data points (stored in t[] and y[],
   containing abscissa and ordinate values respectively), and z[], the
   array of 2nd derivatives at the knots (i.e. internal data points).
   
   On entry, have n+1 >= 2 points in the t, y, z arrays, numbered 0..n.
   The number of knots (i.e. internal data points) is n-1; they are
   numbered 1..n-1 as in Kincaid and Cheney.  In the periodic case, the
   final knot, i.e., (t[n-1],y[n-1]), has the property that y[n-1]=y[0];
   also, y[n]=y[1].  The number of data points supplied by the user was n+1
   in the non-periodic case, and n in the periodic case.  When this
   function is called, n>=1 in the non-periodic case, and n>=2 in the
   periodic case. */

double
#ifdef _HAVE_PROTOS
inline interpolate (int n, double *t, double *y, double *z, double x, 
	     double tension, bool periodic)
#else
inline interpolate (n, t, y, z, x, tension, periodic)
     int n;
     double *t, *y, *z, x;
     double tension;
     bool periodic;
#endif
{
  double diff, updiff, reldiff, relupdiff, h;
  double value;
  int is_ascending = (t[n-1] < t[n]);
  int i = 0, k;

  /* in periodic case, map x to t[0] <= x < t[n] */
  if (periodic && (x - t[0]) * (x - t[n]) > 0.0)
    x -= ((int)(floor( (x - t[0]) / (t[n] - t[0]) )) * (t[n] - t[0]));

  /* do binary search to find interval */
  for (k = n - i; k > 1;)
    {
      if (is_ascending ? x >= t[i + (k>>1)] : x <= t[i + (k>>1)])
	{
	  i = i + (k>>1);
	  k = k - (k>>1);
	}
      else
	k = k>>1;
    }

  /* at this point, x is between t[i] and t[i+1] */
  h = t[i + 1] - t[i];
  diff = x - t[i];
  updiff = t[i+1] - x;
  reldiff = diff / h;
  relupdiff = updiff / h;

  if (tension == 0.0)
  /* evaluate cubic polynomial in nested form */
    {
    value =  y[i] 
      + diff
	* ((y[i + 1] - y[i]) / h - h * (z[i + 1] + z[i] * 2.0) / 6.0
	   + diff * (0.5 * z[i] + diff * (z[i + 1] - z[i]) / (6.0 * h)));
  
    }
  else if (tension > 0.0)
    /* `positive' (really real) tension, use sinh's */
    {
      if (fabs(tension * h) < TRIG_ARG_MIN)
	/* hand-compute (6/y^2)(sinh(xy)/sinh(y) - x) to improve accuracy;
	   here `x' means reldiff or relupdiff and `y' means tension*h */
	value = (y[i] * relupdiff + y[i+1] * reldiff
		 + ((z[i] * h * h / 6.0) 
		    * quotient_sinh_func (relupdiff, tension * h))
		 + ((z[i+1] * h * h / 6.0) 
		    * quotient_sinh_func (reldiff, tension * h)));
      else
	value = (((z[i] * sinh (tension * updiff) 
		   + z[i + 1] * sinh (tension * diff))
		  / (tension * tension * sinh (tension * h)))
		 + (y[i] - z[i] / (tension * tension)) * (updiff / h)
		 + (y[i + 1] - z[i + 1] / (tension * tension)) * (diff / h));
    }
  else
    /* `negative' (really imaginary) tension, use sin's */
    {
     if (fabs(tension * h) < TRIG_ARG_MIN)
	/* hand-compute (6/y^2)(sin(xy)/sin(y) - x) to improve accuracy;
	   here `x' means reldiff or relupdiff and `y' means tension*h */
	value = (y[i] * relupdiff + y[i+1] * reldiff
		 + ((z[i] * h * h / 6.0) 
		    * quotient_sin_func (relupdiff, tension * h))
		 + ((z[i+1] * h * h / 6.0) 
		    * quotient_sin_func (reldiff, tension * h)));
      else
	value = (((z[i] * sin (tension * updiff) 
		   + z[i + 1] * sin (tension * diff))
		  / (tension * tension * sin (tension * h)))
		 + (y[i] - z[i] / (tension * tension)) * (updiff / h)
		 + (y[i + 1] - z[i + 1] / (tension * tension)) * (diff / h));
    }
  
  return value;
}

 
/* is_monotonic() check whether an array of data points, read in by
   read_data(), has monotonic abscissa values. */
bool
#ifdef _HAVE_PROTOS
inline is_monotonic (int n, double *t)
#else
inline is_monotonic (n, t)
     int n;			/* array size n+1, n>=1 */
     double *t;
#endif
{
  bool is_ascending;

  if (t[n-1] < t[n])
    is_ascending = true;
  else if (t[n-1] > t[n])
    is_ascending = false;
  else				/* equality */
    return false;

  while (n>0)
    {
      n--;
      if (is_ascending == true ? t[n] >= t[n+1] : t[n] <= t[n+1])
	return false;
    };
  return true;
}

 

/* Following two functions compute (6/y^2)(sinh(xy)/sinh(y)-x) and
   (6/y^2)(sin(xy)/sin(y)-x), via the first three terms of the appropriate
   power series in y.  They are used when |y|<TRIG_ARG_MIN, to avoid loss
   of significance.  Errors are O(y**6). */
double
#ifdef _HAVE_PROTOS
inline quotient_sinh_func (double x, double y) 
#else
inline quotient_sinh_func (x, y)
     double x, y;
#endif
{
  return ((x*x*x-x) + (x*x*x*x*x/20.0 - x*x*x/6.0 + 7.0*x/60.0)*(y*y)
	  + (x*x*x*x*x*x*x/840.0 - x*x*x*x*x/120.0 + 7.0*x*x*x/360.0
	     -31.0*x/2520.0)*(y*y*y*y));
}

double
#ifdef _HAVE_PROTOS
inline quotient_sin_func (double x, double y) 
#else
inline quotient_sin_func (x, y)
     double x, y;
#endif
{
  return (- (x*x*x-x) + (x*x*x*x*x/20.0 - x*x*x/6.0 + 7.0*x/60.0)*(y*y)
	  - (x*x*x*x*x*x*x/840.0 - x*x*x*x*x/120.0 + 7.0*x*x*x/360.0
	     -31.0*x/2520.0)*(y*y*y*y));
}


void
#ifdef _HAVE_PROTOS
inline non_monotonic_error (void)
#else
inline non_monotonic_error ()
#endif
{
  fprintf (stderr, "%s: error: abscissa values not monotonic\n",
	   progname);
  exit (EXIT_FAILURE);
}

 


/* Following four functions compute (6/x^2)(1-x/sinh(x)),
   (3/x^2)(x/tanh(x)-1), (6/x^2)(1-x/sin(x)), and (3/x^2)(x/tan(x)-1) via
   the first three terms of the appropriate power series.  They are used
   when |x|<TRIG_ARG_MIN, to avoid loss of significance.  Errors are
   O(x**6). */
double
#ifdef _HAVE_PROTOS
inline sinh_func (double x) 
#else
inline sinh_func (x)
     double x;
#endif
{
  /* use 1-(7/60)x**2+(31/2520)x**4 */
  return 1.0 - (7.0/60.0)*x*x + (31.0/2520.0)*x*x*x*x;
}

double
#ifdef _HAVE_PROTOS
inline tanh_func (double x) 
#else
inline tanh_func (x)
     double x;
#endif
{
  /* use 1-(1/15)x**2+(2/315)x**4 */
  return 1.0 - (1.0/15.0)*x*x + (2.0/315.0)*x*x*x*x;
}

double
#ifdef _HAVE_PROTOS
inline sin_func (double x) 
#else
inline sin_func (x)
     double x;
#endif
{
  /* use -1-(7/60)x**2-(31/2520)x**4 */
  return -1.0 - (7.0/60.0)*x*x - (31.0/2520.0)*x*x*x*x;
}

double
#ifdef _HAVE_PROTOS
inline tan_func (double x) 
#else
inline tan_func (x)
     double x;
#endif
{
  /* use -1-(1/15)x**2-(2/315)x**4 */
  return -1.0 - (1.0/15.0)*x*x - (2.0/315.0)*x*x*x*x;
}

 
