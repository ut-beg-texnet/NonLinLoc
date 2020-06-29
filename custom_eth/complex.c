#include <math.h>

typedef struct FCOMPLEX {float r,i;} fcomplex;
typedef struct DCOMPLEX {double r,i;} dcomplex;



dcomplex Cadd(a,b)
dcomplex a,b;
{
	dcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

dcomplex Csub(a,b)
dcomplex a,b;
{
	dcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}


dcomplex Cmul(a,b)
dcomplex a,b;
{
	dcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

dcomplex Complex(re,im)
float im,re;
{
	dcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

dcomplex Conjg(z)
dcomplex z;
{
	dcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

dcomplex Cdiv(a,b)
dcomplex a,b;
{
	dcomplex c;
	double r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

double Cabs(z)
dcomplex z;
{
	double x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

dcomplex Csqrt(z)
dcomplex z;
{
	dcomplex c;
	double x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=fabs(z.r);
		y=fabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

dcomplex RCmul(x,a)
dcomplex a;
float x;
{
	dcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}
