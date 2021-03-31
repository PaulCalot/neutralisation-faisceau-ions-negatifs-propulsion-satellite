#include <stdio.h>
#include <math.h>

int build_=1;
void main (int argc, char * argv[])
{
    double Vmol_ra, V1mol_ra, V2mol_ra, Vmor_rb, V1mor_rb, V2mor_rb;
    double morseA, morseB, morseLam, morseMu, moliA, moliC;
    double Vr, dVr, d2Vr;
    double a, b, c, s;
    double ra, rb, r;
    void MorseVr(double, double, double, double *, double *);
    void MoliVr(double, double, double, double *, double *);
    

    if (argc<7)
    {
	printf("Usage: coeff3 ra rb morseA morseLam moliA moliC\n");
	exit(-1);
    }
    
    ra = atof(argv[1]);
    rb = atof(argv[2]);
    morseA = atof(argv[3]);
    morseLam = atof(argv[4]);
    moliA = atof(argv[5]);
    moliC = atof(argv[6]);
    
    MorseVr(rb, morseA, morseLam, &Vmor_rb, &V1mor_rb);
    MoliVr(ra, moliA, moliC, &Vmol_ra, &V1mol_ra);

    a = 1.0/(ra-rb)*log(V1mol_ra/V1mor_rb);
    b = log(V1mol_ra/a) - a*ra;
    c = Vmor_rb - exp(a*rb+b);
    s = exp(a*ra+b) + c - Vmol_ra;

    printf("#define MS_RA\t    %.10le\n", ra);
    printf("#define MS_RB\t    %.10le\n", rb);
    printf("#define MS_A\t    %.10le\n", a);
    printf("#define MS_B\t    %.10le\n", b);
    printf("#define MS_C\t    %.10le\n", c);
    printf("#define MS_S\t    %.10le\n", s);
}

void MorseVr (double r, double A, double lam, 
    double * Vr0, double * Vr1)
{
    *Vr0 = A*exp(-lam*r);
    *Vr1 = (*Vr0)*(-lam);    
}

static const double m_c[3] = {0.35, 0.55, 0.1};
static const double m_d[3] = {0.3, 1.2, 6.2};
static const double m_cd[3] = {0.105, 0.66, 0.62};
static const double m_cdd[3] = {0.0315, 0.792, 3.844};
static double rij_am, m_efac, m_v1_s1, m_v1_s2, m_v1_s3, mc_efac;
void MoliVr (double rij, double A, double C, double * Vr0, double * Vr1)
{
    int i = 0;
    
    *Vr0 = *Vr1 = m_v1_s1 = m_v1_s2 = m_v1_s3 =  0.0;
    rij_am = rij/A;	/* dimensionless rij for the Moliere function */
    
    for (i = 0; i < 3; i++)
    {
	m_efac = exp(-m_d[i]*rij_am);
	mc_efac = m_c[i]*m_efac;
	*Vr0 += mc_efac;
	m_v1_s1 += mc_efac;
	m_v1_s2 += m_cd[i]*m_efac;
    }
    *Vr0 *= C/rij_am;
    m_v1_s1 /= rij;
    m_v1_s2 /= A;
    *Vr1 = m_v1_s1 + m_v1_s2;
    *Vr1 *= -C/rij_am;
}
