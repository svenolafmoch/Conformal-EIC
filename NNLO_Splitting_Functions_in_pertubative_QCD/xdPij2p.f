*
* ..File: xdPij2p.f
*
*                         __
* ..The parametrized NNLO MS (M version, as defined in hep-ph/9803439)
*    flavour-singlet splitting functions  DeltaP^(2)  for the evolution
*    of helicity-difference (polarized) parton densities at mu_r = mu_f.
*    The expansion parameter is alpha_s/(4 pi).
*
* ..The distributions (in the mathematical sense) are given as in eq.
*    (B.27) of Floratos, Kounnas, Lacaze: Nucl. Phys. B192 (1981) 417.
*    The name-endings A, B, and C of the functions below correspond to
*    the kernel superscripts [2], [3], and [1] in that equation.
*
* ..The relative accuracy of these parametrisations, as well as of the
*    resulting convolutions, is better than one part in a thousand.
*
* ..The coefficients of ln^4(x) are exact, as are those of 1/(1-x)_+
*    in DeltaP_gg^(2) up to a truncation of irrational coefficients.
*    The other terms at x < 1 have fitted to the exact results for x
*    between 10^-6 and 1 - 10^-6.  The coefficients of delta(1-x) of
*    DeltaP_gg^(2) have been slightly adjusted using low moments.
*
* ..Reference: S. Moch, J. Vermaseren and A. Vogt,
*                The Three-Loop Splitting Functions in QCD: 
*                The Helicity-Dependent Case 
*              September 2014
*
* =====================================================================
*
*
* ..The pure-singlet splitting function DeltaP_ps^(2). 
*    The quark-quark splitting function DeltaP_qq^(2) is obtained by 
*    adding the unpolarized quantity P_ns^(2)- given in hep-ph/0403192.

       FUNCTION DP2PSPA (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf
*
       L  = log (x)
       L1 = log (1.d0-x)
*
       P2ps1 = - 344./27.d0 * L**4 - (90.9198 + 81.50* x)* L**3 
     ,         - (368.6 - 349.9* x)* L*L - (739.0 - 232.57* L1)* L 
     ,         - 1362.6 + 1617.4 * x - 674.8 * x*x + 167.41 * x**3
     ,         - 204.76 * L1 - 12.61 * L1*L1 - 6.541 * L1**3
       P2ps2 = (1.1741 - 0.8253* x)* L**3  + (13.287 + 10.657* x)* L*L 
     ,         + 45.482 * L + 49.13 - 30.77 * x - 4.307 * x*x 
     ,         - 0.5094 *x**3 + 9.517 * L1 + 1.7805 * L1*L1
*
       DP2PSPA = (1.-x) * nf * ( P2ps1 + nf * P2ps2 ) 
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon->quark splitting function DeltaP_qg^(2).
*
       FUNCTION DP2QGPA (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf 
*
       L  = log (x)
       L1 = log (1.d0-x)
*
       P2qg1 = - 151./3.d0 * L**4 - (385.64 + 73.30* x)* L**3
     ,         - (894.8 - 1145.3* x)* L*L - (1461.2 - 825.4* L1)* L
     ,         - 2972.4 + 4672.* x - 1221.6 * x*x - 18.0 * x**3
     ,         + 278.32* L1 - 90.26* L1*L1 - 5.30* L1**3 + 3.784*L1**4
       P2qg2 =   16./9.d0 * L**4 + (30.739  + 10.186* x) * L**3 
     ,         + (196.96 + 179.1* x)* L*L + (526.3  - 47.30* L1)* L 
     ,         + 499.65 - 432.18 * x - 141.63 * x*x - 11.34 * x**3 
     ,         - 6.256 * L1 + 7.32 * L1*L1 + 0.7374 * L1**3
* 
       DP2QGPA = nf * ( P2qg1 + nf * P2qg2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The quark->gluon splitting functions DeltaP_gq^(2). P2gq2 is exact.
*
       FUNCTION DP2GQPA (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf 
*
       L  = log (x)
       L1 = log (1.d0-x)
*
       P2gq0 =   11512./81.d0* L**4 + (888.003  + 175.1* x)* L**3
     ,         + (2140. - 850.7* x)* L*L + (4046.6 - 1424.8* L1)* L
     ,         + 6159. - 3825.9 * x + 1942.* x*x - 742.1 * x**3
     ,         + 1843.7* L1 + 451.55* L1*L1 + 59.3* L1**3 + 5.143* L1**4
       P2gq1 = - 128./27.d0 * L**4 - (39.3872 + 30.023*x)* L**3
     ,         - (202.46 + 126.53* x)* L*L - (308.98 + 16.18* L1)* L
     ,         - 301.07 - 296.0 * x + 406.13 * x*x - 101.62 * x**3
     ,         - 171.78* L1 - 47.86 * L1*L1 - 4.963 * L1**3
       P2gq2 =   16./27.d0 * ( - 12. + 10.* x + ( 8.+ 2.*x)* L1 
     ,         + (6.- 3.*x)* L1*L1 )
*
       DP2GQPA = P2gq0 + nf * (P2gq1 + nf * P2gq2) 
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon-gluon splitting function DeltaP_gg^(2), regular piece.
*
       FUNCTION DP2GGPA (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER NF
*
       L  = log (x)
       L1 = log (1.d0-x)
*
       P2ggA0 =   504.d0 * L**4 + (3777.5  + 1167.* x)* L**3
     ,          + (10902. - 863.* x)* L*L + (23091. - 12292.* L1)* L
     ,          + 30988. - 39925.* x + 13447.* x*x - 4576.* x**3 
     ,          - 13247.* (1.-x)*L1 + 3801.* L1
       P2ggA1 = - 766./27.d0 * L**4 - (357.798 - 131.* x)* L**3
     ,          - (1877.2 - 613.1* x)* L*L - (3524. + 7932.* L1)* L
     ,          - 1173.5 + 2648.6 * x - 2160.8 * x*x + 1251.7 * x**3
     ,          - 6746.* (1.-x)*L1 - 295.7* L1
       P2ggA2 = - 1.1809 * L**3 - (6.679 - 15.764* x)* L*L 
     ,          - (13.29 + 16.944* L1) * L - 16.606 + 32.905 * x 
     ,          - 18.30 * x*x + 2.637 * x**3 - 0.210 * L1
*
       DP2GGPA = P2ggA0 + nf * ( P2ggA1 + nf * P2ggA2 )
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon-gluon splitting function DeltaP_gg^(2), singular piece.
*
       FUNCTION DP2GGPB (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf
*
       DP2GGPB = ( 2643.521 - nf * (412.172 + nf * 16./9.D0 ) ) 
     1           / (1.d0-x)
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The gluon-gluon splitting function DeltaP_gg^(2), `local' piece.  
*
       FUNCTION DP2GGPC (x, nf)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf
*
       L1 = log (1.d0-x)
*
       DP2GGPC =             2643.521 * L1 + 4425.448 + 2.314   
     ,           - nf *    (  412.172 * L1 +  528.720 - 0.184 ) 
     ,           - nf*nf * ( 16./9.D0 * L1 -   6.4630 + 0.0023 )
*
       RETURN
       END
*
* =================================================================av==
