/* 
 * md series 2 (c) 1999 Cameron Abrams
 *
 * sicf_params.h
 * 
 * Parameters for the Tersoff/Brenner SiCF interatomic potential as given in:
 * 
 * (1) Beardmore96 (PhilMagA      74(6) p1439 (1996))
 * (2) Brenner90   (PhysRevB     42(15) p9458 (1990))
 * (3) Tanaka99    (JVacSciTechA accepted (1999))
 * (4) Abrams99c   (JApplPhys    86(11) p5938 (1999))
 *
 * Note:  Unless otherwise noted, all parameter values in this
 * file are from Parameter set B.  Those that end in __A are from
 * parameter set A (the older set).
 *
 */

/* Tersoff, Moliere Repulsive, and Moliere-to-Tersoff Spline 
 * parameter sets.  The Tersoff-form parameters are computed
 * from Brenner-form parameters using m2t.c (c) 1999 C. Abrams.
 * Parameters for the Moliere-to-Tersoff Spline computed using
 * coeff3.c (c) 1999 C. Abrams. */

/* Explanation of two-body potential parameters:
 * 
 * The 4-parameter Tersoff-form 2body potential is:
 * V(r) = A*exp(-lamda*r) - B*exp(-mu*r)
 * 
 * The 4-parameter Morse/Brenner-form 2body potential is:
 * Vmb(r) = (D/(S-1))*exp(-sqrt(2*S)*Beta*(r-Re)) - 
 *	    (D*S/(S-1))*exp(-sqrt(2/S)*Beta*(r-Re)) 
 * 
 * Note:  The Tersoff-form is adopted to describe all bonds in
 * this implementation.  As the Tersoff and Morse/Brenner forms
 * are functionally identical wrt interatomic separation "r", 
 * we are free to translate the 4 M/B form paramters into the
 * 4 Tersoff-form parameters.  This is what I have done so that
 * Tersoff-form PEF is used consistently for all types of bonds.
 * This follows the method of Beardmore96.
 *
 * The 3-parameter Moliere repulsive potential is :
 * Vm(r) = C/r*sum(i=1,4)[c_i*exp(-d_i*r/a)] + s, 
 *   where C, a and s are the parameters.  (C is actually the
 *   nuclear-nuclear Coulomb factor, (ZiZje^2)/(4*pi*ep_0).)
 * 
 * The 3-parameter Moliere-to-Tersoff spline potential is:
 * Vs(r) = c + exp(a*r + b), 
 *   where c, a, and b are the parameters.
 */

/* Si-Si 2body parameter values */
/* Tersoff */
#define A_SiSi	    1830.800		/* eV */
#define B_SiSi	    471.1800		/* eV */
#define LAMBDA_SiSi 2.479900		/* 1/Angstroms */
#define MU_SiSi	    1.732200		/* 1/Angstroms */
/* Moliere-to-Tersoff Spline, including the Moliere "s" parameter */
#define MS_RASiSi   0.286968		/* Angstroms */
#define MS_RBSiSi   0.652200		/* Angstroms */
#define MS_ASiSi    -7.1553762043e+00	/* 1/Angstroms */
#define MS_BSiSi    9.5022081352e+00
#define MS_CSiSi    2.3736156216e+02	/* eV */
#define MS_SSiSi    2.9679293202e+02	/* eV */
/* Moliere parameters "a" and "C" */
#define M_ASiSi	    0.10164		/* Angstroms */
#define M_CSiSi	    27769.25601141282959464777	    /* eV */

#define RSISI_CU1	2.7		/* R(1), Angstroms */
#define RSISI_CU2	3.0		/* R(2), Angstroms */
#define RSISI_CU2SQ	9.0		/* R(2)*R(2), sqAngstroms */
#define PI_RSISI_CU12	10.47197551	/* pi/(R(2)-R(1)), Angstroms */
#define RSISI_E		2.35		/* R(e), Angstroms */

/* The CC Tersoff parameters are computed from the Morse/Brenner-form
 * parameters reported in Brenner90.  These paramter values belong both
 * to set A and set B.
 * 
 * Morse/Brenner-form parameters: From Table II in above reference
 * R(e)(CC) = 1.39 A
 * D(CC)    = 6.0 eV
 * Beta(CC) = 2.1   1/A
 * S(CC)    = 1.22
 * 
 * Tersoff-form parameters as determined from M/B-form parameters:
 * A = D/(S-1)*exp(sqrt(2*S)*R(e)*Beta)	    = 2605.8416 eV
 * B = D*S/(S-1)*exp(sqrt(2/S)*Beta*R(e))   = 1397.0730 eV
 * lambda = sqrt(2*S)*Beta		    = 3.2803 1/A
 * mu = sqrt(2/S)*Beta			    = 2.6888 1/A
 */
#define A_CC	    2605.8416	    /* eV */
#define B_CC	    1397.0730	    /* eV */
#define LAMBDA_CC   3.2803	    /* 1/Angstroms */
#define MU_CC	    2.6888	    /* 1/Angstroms */
#define MS_RACC	    2.8696800000e-01
#define MS_RBCC	    6.5220000000e-01	/* Beardmore uses 0.369580 */
#define MS_ACC	    -2.8628849403e+00
#define MS_BCC	    7.7293779159e+00
#define MS_CCC	    -4.4727801734e+01
#define MS_SCC	    5.4423766652e+02
#define M_ACC	    0.13481			    /* Angstroms */
#define M_CCC	    3845.50359320525183591721	    /* eV */

#define RCC_CU1		1.7		/* Angstroms */
#define RCC_CU2		2.0		/* Angstroms */
#define RCC_CU2SQ	4.0		/* sqAngstroms */
#define PI_RCC_CU12	10.47197551	/* Angstroms */
#define RCC_E		1.39		/* Angstroms */


/* The FF Tersoff Parameters are computed from the Morse/Brenner-form
 * parameters reported in Tanaka99.
 * 
 * Brenner-form parameters:
 * R(e)(FF) = 1.4119 A
 * D(FF)    = 1.6020 eV
 * Beta(FF) = 3.026  1/A
 * S(FF)    = 1.176
 * 
 * Tersoff-form parameters as determined from M/B-form parameters:
 * A = D/(S-1)*exp(sqrt(2*S)*R(e)*Beta)	    = 6379.1500 eV
 * B = D*S/(S-1)*exp(sqrt(2/S)*Beta*R(e))   = 2813.8185 eV
 * lambda = sqrt(2*S)*Beta		    = 4.6407 1/A
 * mu = sqrt(2/S)*Beta			    = 3.9462 1/A
 */
#define A_FF	    6379.1500	    /* eV */
#define B_FF	    2813.8185	    /* eV */
#define LAMBDA_FF   4.6407 	    /* 1/Angstroms */
#define MU_FF	    3.9462	    /* 1/Angstroms */
#define MS_RAFF	    2.8696800000e-01
#define MS_RBFF	    6.5220000000e-01
#define MS_AFF      -3.8113458035e+00
#define MS_BFF      8.4167681926e+00
#define MS_CFF      -6.7291591284e+01
#define MS_SFF      6.4252174315e+02
#define M_AFF	    0.11777	    /* Angstroms */
#define M_CFF	    9904.28601511420565509043	    /* eV */

/* parameter set A */
#define A_FF__A		16451.97	    /* eV */ 
#define B_FF__A		146.8149	    /* eV */ 
#define LAMBDA_FF__A    6.8149		    /* 1/Angstroms */ 
#define MU_FF__A        2.8568		    /* 1/Angstroms */ 
#define MS_AFF__A       -4.0477892098e+00
#define MS_BFF__A       8.4244313820e+00
#define MS_CFF__A       -1.3204669470e+02
#define MS_SFF__A       4.8928309004e+02           
#define RFF_E__A	1.4119		/* Angstroms */

#define RFF_CU1		1.7		/* Angstroms */
#define RFF_CU2		2.0		/* Angstroms */
#define RFF_CU2SQ	4.0		/* sqAngstroms */
#define PI_RFF_CU12	10.47197551	/* Angstroms */
#define RFF_E		1.4119		/* Angstroms */

/* SiC parameters taken directly from Beardmore96. */
#define A_SiC	    1597.311	    /* eV */
#define B_SiC	    395.1451	    /* eV */
#define LAMBDA_SiC  2.983900	    /* 1/Angstroms */
#define MU_SiC	    1.972050	    /* 1/Angstroms */
#define MS_RASiC    0.286968	    /* Angstroms */
#define MS_RBSiC    0.652200	    /* Angstroms */
#define MS_ASiC	    -5.90430	    /* 1/Angstroms */
#define MS_BSiC	    8.59831
#define MS_CSiC	    112.845	    /* eV */
#define MS_SSiC	    292.680	    /* eV */
#define M_ASiC	    0.11533	    /* Angstroms */
#define M_CSiC	    10488.41404664874707361484     /* eV */

#define RSIC_CU1	2.204541	/* Angstroms */
#define RSIC_CU2	2.509980	/* Angstroms */
#define RSIC_CU2SQ	6.299999	/* sqAngstroms */
#define PI_RSIC_CU12	10.28549941	/* Angstroms */
#define RSIC_E		1.85		/* Angstroms */

/* The CF Tersoff Parameters are computed from the M/B-form
 * parameters reported in Tanaka99.
 * 
 * Morse/Brenner-form parameters:
 * R(e)(CF) = 1.2667 A
 * D(CF)    = 5.5642 eV
 * Beta(CF) = 2.040 1/A
 * S(CF)    = 1.1978	= sqrt[S(FF)*S(CC)] = sqrt(1.176*1.22)
 * 
 * Tersoff-form parameters as determined from Morse/Brenner-form parameters:
 * A = D/(S-1)*exp(sqrt(2*S)*R(e)*Beta)	    = 1535.1781 eV
 * B = D*S/(S-1)*exp(sqrt(2/S)*Beta*R(e))   = 949.9585 eV
 * lambda = sqrt(2*S)*Beta		    = 3.1575 1/A
 * mu = sqrt(2/S)*Beta			    = 2.6360 1/A
 */
#define A_CF	    1535.1781	    /* eV */
#define B_CF	    949.9585	    /* eV */
#define LAMBDA_CF   3.1575	    /* 1/Angstroms */
#define MU_CF	    2.6360	    /* 1/Angstroms */
#define MS_RACF	    2.8696800000e-01
#define MS_RBCF	    6.5220000000e-01
#define MS_ACF      -5.1521258108e+00
#define MS_BCF      8.1476678001e+00
#define MS_CCF      7.5802194262e+01
#define MS_SCF      2.8973644378e+02
#define M_ACF	    0.12557	    /* Angstroms */
#define M_CCF	    6192.70931751214462053038	    /* eV */

/* parameter set A */
#define A_CF__A		909.2022	    /* eV */ 
#define B_CF__A		219.7799	    /* eV */ 
#define LAMBDA_CF__A    3.7128		    /* 1/Angstroms */ 
#define MU_CF__A        2.1463		    /* 1/Angstroms */ 
#define MS_ACF__A	-7.1344140384
#define MS_BCF__A	8.3910003739
#define MS_CCF__A	3.8716335228e+01
#define MS_SCF__A	3.3777008484e+01
#define RCF_E__A	1.2718		/* Angstroms */

#define RCF_CU1		1.7		/* Angstroms */
#define RCF_CU2		2.0		/* Angstroms */
#define RCF_CU2SQ	4.0		/* sqAngstroms */
#define PI_RCF_CU12	10.47197551	/* Angstroms */
#define RCF_E		1.2667		/* Angstroms */

/* SiF parameters as reported in Abrams99c. */

#define A_SiF	    37412.28	    /* eV */
#define B_SiF	    925.846	    /* eV */
#define LAMBDA_SiF  5.4875	    /* 1/Angstroms */
#define MU_SiF	    2.7437	    /* 1/Angstroms */
#define MS_RASiF    3.0000000000e-01
#define MS_RBSiF    8.0000000000e-01
#define MS_ASiF	    -2.1325772404e+00
#define MS_BSiF	    8.7909576222e+00
#define MS_CSiF	    -7.2985931698e+02
#define MS_SSiF	    1.6885525709e+03
#define M_ASiF	    0.10897	    /* Angstroms */
#define M_CSiF	    16650.85058272919152060200	    /* eV */

#define RSIF_CU1	1.83922		/* Angstroms */
#define RSIF_CU2	2.13922		/* Angstroms */
#define RSIF_CU2SQ	4.576262208	/* sqAngstroms */
#define PI_RSIF_CU12	10.47197551	/* Angstroms */
#define RSIF_E		1.6008		/* Angstroms */

/* Moliere repulsive parameters for (Si,C,F)-Ar */
#define M_ASiAr	    0.09734 /* Angstroms */
#define M_CSiAr	    37280.52574481199917813848 /* eV */
#define MS_SSiAr    0.00000
#define MS_RASiAr   1.e10     /* ensures Moliere is always used */
#define MS_RBSiAr   MS_RASiAr /* ensures spline is never used */
#define RSIAr_CU2SQ 25.0      /* ensures a 5-Angstrom cutoff of Ar-Si */

#define M_ACAr	    0.10950 /* Angstroms */
#define M_CCAr	    14203.07778995433789954337 /* eV */
#define MS_SCAr     0.00000
#define MS_RACAr    1.e10
#define MS_RBCAr    MS_RACAr
#define RCAr_CU2SQ  16.0      /* ensures a 4-Angstrom cutoff of Ar-C */

#define M_AFAr	    0.10388 /* Angstroms */
#define M_CFAr	    22457.21531574894108586830 /* eV */
#define MS_SFAr     0.00000
#define MS_RAFAr    1.e10
#define MS_RBFAr    MS_RACAr
#define RFAr_CU2SQ  25.0      /* ensures a 5-Angstrom cutoff of Ar-F */

/* Three-body parameters */

#define eta_Cbrenner	    1.0
#define eta_Ctersoff	    0.72752
#define eta_Sitersoff	    0.78734
#define eta_Simurty	    1.0
#define eta_Ftanaka	    1.0		/* implicit in Tanaka99 */
#define delta_Cbrenner	    0.5
#define delta_Ctersoff	    0.687276
#define delta_Sitersoff	    0.635050
#define delta_Simurty	    0.80469
#define delta_Ftanaka	    0.5	    	/* in Tanaka99 */
#define alpha_Cbrenner	    6.0		/* in Tanaka99 */
#define alpha_Cbrenner__A   3.0
#define alpha_Sitersoff	    5.197495	/* 1/Angstrom */
#define alpha_Simurty	    4.0		/* 1/Angstrom */
#define alpha_Ftanaka	    6.0		/* in Tanaka99 */
#define alpha_Ftanaka__A    3.0
#define beta_Cbrenner	    1
#define beta_Sitersoff	    3
#define beta_Simurty	    3
#define beta_Ftanaka	    1		/* implicit in Tanaka99 */

/* Parameters for the g(cosTheta) functions */
#define a_CbrennerI	0.011304
#define a_Cbrenner	0.00020813
#define a_Ctersoff	1.5724e-7
#define c_CbrennerI	19.0
#define c_Cbrenner	330.0
#define c2_CbrennerI	361.0		/* c^2 for the C-brennerI g() */
#define c2_Cbrenner	108900.0	/* c^2 for the C-brenner g() */
#define ac2_CbrennerI	4.080744	/* a*c^2 for the C-brenner g() */
#define ac2_Cbrenner	22.665357	/* a*c^2 for the C-brenner g() */
#define c_Ctersoff	3.8049e-4
#define c2_Ctersoff	1.4477264e-7	/* c^2 for the C-tersoff g() */
#define ac2_Ctersoff	2.276404991e-14	/* a*c^2 for the C-tersoff g() */
#define c_Sitersoff	0.0
#define d_CbrennerI	2.5
#define d2_CbrennerI	6.25		/* d^2 for the C-brenner g() */
#define c2d2_CbrennerI	57.76		/* c^2/d^2 for the C-brenner g() */
#define ac2d2_CbrennerI	0.65291904	/* a*c^2/d^2 for the C-brenner g() */
#define d_Cbrenner	3.5
#define d2_Cbrenner	12.25		/* d^2 for the C-brenner g() */
#define c2d2_Cbrenner	8889.795918	/* c^2/d^2 for the C-brenner g() */
#define ac2d2_Cbrenner	1.850233224	/* a*c^2/d^2 for the C-brenner g() */
#define d_Ctersoff	4.384
#define d2_Ctersoff	19.219456	/* d^2 for the C-tersoff g() */
#define c2d2_Ctersoff	7.532608623e-9	/* c^2/d^2 for the C-tersoff g() */
#define ac2d2_Ctersoff	1.18442738e-15	/* a*c^2/d^2 for the C-tersoff g() */
#define d_Sitersoff	0.16
#define c_Simurty	0.0216
#define d_Simurty	0.27
#define h_Simurty	-0.470
#define h_Cbrenner	-1.0
#define h_Ctersoff	-0.57058
#define h_Sitersoff	-0.59826

#define a_Ftanaka	0.0 /* not used */
#define c_Ftanaka__A	14.0  /* so G returns a constant */
#define c_Ftanaka	4.0  /* so G returns a constant */
#define d_Ftanaka	0.0 /* not used */
#define h_Ftanaka	0.0 /* not used */

/* Precomputed factors involving parameters of the g(cosTheta) functions */
#define p2_CbrennerI	ac2_CbrennerI
#define p1_CbrennerI	0.66422304
#define p2_Cbrenner	ac2_Cbrenner
#define p1_Cbrenner	1.850441354
#define p2_Ctersoff	ac2_Ctersoff
#define p1_Ctersoff	1.5724000118442738e-7
