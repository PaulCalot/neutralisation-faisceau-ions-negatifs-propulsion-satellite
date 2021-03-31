/* converts the 4 brenner/morse-form paramters to 
 * tersoff/beardmore-form parameters:
 *
 *  Brenner/Morse		Tersoff
 *  R_e  bond length,		A
 *  D_e  bond strength,		B
 *  Beta vib.,			lam
 *  S    something else		mu
 */
#include <stdio.h>
#include <math.h>

int main (int argc,char * argv[])
{
	double r,d,b,s;
	double a,B,l,m;

	r=atof(argv[1]);
	d=atof(argv[2]);
	b=atof(argv[3]);
	s=atof(argv[4]);

	a=d/(s-1)*exp(sqrt(2*s)*r*b);
	B=d*s/(s-1)*exp(sqrt(2./s)*r*b);
	l=sqrt(2*s)*b;
	m=sqrt(2./s)*b;

	printf("%.8lf %.8lf %.8lf %.8lf\n",
		a,B,l,m);
}
