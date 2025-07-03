*
* ..File: xptsdiff.f 
*
*
* ..The differences between the time-like and space-like non-singlet 
*    splitting functions, respectively governing the evolution of the
*    fragmentation functions and the parton distributions, at second 
*    and third order in alpha_s/(4 pi) for the scale mu_f = mu_r.
* 
* ..The code uses the package of Gehrmann and Remiddi for the harmonic
*    polylogarithms published in hep-ph/0107173 = CPC 141 (2001) 296.
*
* ..Reference: A. Mitov, S. Moch and A. Vogt,  hep-ph/0604053
*
* ..These functions need to be combined with the space-like results of
*              S. Moch, J. Vermaseren and A. Vogt, 
*              hep-ph/0403192 = Nucl. Phys. B688 (2004) 101
*
* =====================================================================
*
*
* ..The second-order (NLO) function P_NS^T(1) - P_NS^S(1) for all
*   non-singlet combinations.
*
       FUNCTION DTSP1NS (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       INTEGER NF
*
* ..Colour factor and abbreviations
*
       CF = 4.D0/3.D0
       DM = 1.D0/(1.D0-X)
       LX = LOG(X)
*
* ..The splitting-function difference in terms of logarithms 
*
       DTSP1NS = LX*LOG(1.D0-X) * ( 32.D0* DM - 16.D0 - 16.D0* X )
     ,         + LX**2 * (  - 16.D0* DM + 12.D0 + 12.D0* X )
     ,         + LX* ( 24.D0* DM - 20.D0 - 4.D0* X );
       DTSP1NS = CF**2 * DTSP1NS
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The third-order (NNLO) function P_NS+^T(2) - P_NS+^S(2) for the 
*   flavour differences of quark-antiquark sums.  
*
       FUNCTION DTSP2NSP (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4 
       INTEGER NF, N1, N2, NW, I1, I2, I3, N
       PARAMETER ( N1 = -1, N2 = 1, NW = 4 ) 
       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
*
* ...Colour factors
*
       CF  = 4./3.D0
       CA  = 3.D0
       CAM2CF = CA - 2.D0 * CF
*
* ...Some abbreviations
*
       PQQ0P = 2.D0/(1.D0-X) - 1.D0 - X
       PQQ0M = 2.D0/(1.D0+X) - 1.D0 + X
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2) 
*
* ...The splitting-function difference in terms of the harmonic polylogs
*
      ppdiff =
     &  + cf**2*cam2cf * (  - 88.D0/9.D0*Hr1(0) + 88.D0/9.D0*Hr1(0)*x
     &     + 3.D0*Hr1(0)*z3 - 3.D0*Hr1(0)*z3*x - 4.D0*Hr1(0)*z2*x - 50.D
     &    0/9.D0*Hr2(0,0) + 184.D0/9.D0*Hr2(0,0)*x - 4.D0*Hr2(0,0)*z2*x
     &     - 28.D0/3.D0*Hr2(0,1) + 28.D0/3.D0*Hr2(0,1)*x - 28.D0/3.D0*
     &    Hr2(1,0) + 28.D0/3.D0*Hr2(1,0)*x - 8.D0*Hr3(-1,0,0) - 8.D0*
     &    Hr3(-1,0,0)*x - 4.D0*Hr3(0,-1,0) - 4.D0*Hr3(0,-1,0)*x + 11.D0/
     &    2.D0*Hr3(0,0,0) + 35.D0/2.D0*Hr3(0,0,0)*x + 4.D0*Hr4(0,-1,0,0
     &    ) - 4.D0*Hr4(0,-1,0,0)*x + 4.D0*Hr4(0,0,-1,0) - 4.D0*Hr4(0,0,
     &    -1,0)*x + 8.D0*Hr4(0,0,0,0)*x )
      ppdiff = ppdiff + cf**2*cam2cf*pqq0m * (  - 7.D0*Hr1(0)*z3 + 3.D
     &    0/2.D0*Hr1(0)*z2 + 8.D0*Hr2(-1,0)*z2 + 8.D0*Hr2(0,-1)*z2 - 8.D
     &    0*Hr2(0,0)*z2 + 6.D0*Hr3(-1,0,0) + 3.D0*Hr3(0,-1,0) - 9.D0/2.D
     &    0*Hr3(0,0,0) + 16.D0*Hr4(-1,-1,0,0) + 8.D0*Hr4(-1,0,-1,0) - 
     &    18.D0*Hr4(-1,0,0,0) - 8.D0*Hr4(-1,0,0,1) - 4.D0*Hr4(-1,0,1,0)
     &     + 8.D0*Hr4(0,-1,-1,0) - 14.D0*Hr4(0,-1,0,0) - 4.D0*Hr4(0,-1,
     &    0,1) - 8.D0*Hr4(0,0,-1,0) + 8.D0*Hr4(0,0,0,0) + 6.D0*Hr4(0,0,
     &    0,1) + 2.D0*Hr4(0,0,1,0) )
      ppdiff = ppdiff + cf**2*cam2cf*pqq0p * ( 151.D0/24.D0*Hr1(0) + 
     &    Hr1(0)*z3 + 13.D0/6.D0*Hr1(0)*z2 - 169.D0/18.D0*Hr2(0,0) + 8.D
     &    0*Hr2(0,0)*z2 - 134.D0/9.D0*Hr2(0,1) + 4.D0*Hr2(0,1)*z2 - 134.
     &    D0/9.D0*Hr2(1,0) + 4.D0*Hr2(1,0)*z2 - 13.D0/2.D0*Hr3(0,0,0)
     &     - 22.D0/3.D0*Hr3(0,0,1) - 22.D0/3.D0*Hr3(0,1,0) - 22.D0/3.D0
     &    *Hr3(1,0,0) - 8.D0*Hr4(0,0,0,0) - 6.D0*Hr4(0,0,0,1) - 2.D0*
     &    Hr4(0,0,1,0) - 2.D0*Hr4(0,1,0,0) - 6.D0*Hr4(1,0,0,0) )
      ppdiff = ppdiff + cf**3 * (  - 325.D0/18.D0*Hr1(0) + 325.D0/18.D0
     &    *Hr1(0)*x + 3.D0*Hr1(0)*z2 - 5.D0*Hr1(0)*z2*x - 173.D0/18.D0*
     &    Hr2(0,0) + 691.D0/18.D0*Hr2(0,0)*x - 4.D0*Hr2(0,0)*z2 - 4.D0*
     &    Hr2(0,0)*z2*x - 50.D0/3.D0*Hr2(0,1) + 50.D0/3.D0*Hr2(0,1)*x
     &     - 50.D0/3.D0*Hr2(1,0) + 50.D0/3.D0*Hr2(1,0)*x + 25.D0/2.D0*
     &    Hr3(0,0,0) + 25.D0/2.D0*Hr3(0,0,0)*x + 2.D0*Hr3(0,0,1) + 2.D0
     &    *Hr3(0,0,1)*x + Hr3(0,1,0) + Hr3(0,1,0)*x )
      ppdiff = ppdiff + cf**3*pqq0p * ( 311.D0/24.D0*Hr1(0) + 4.D0/3.D0
     &    *Hr1(0)*z2 - 169.D0/9.D0*Hr2(0,0) + 8.D0*Hr2(0,0)*z2 - 268.D0/
     &    9.D0*Hr2(0,1) + 8.D0*Hr2(0,1)*z2 - 268.D0/9.D0*Hr2(1,0) + 8.D0
     &    *Hr2(1,0)*z2 - 22.D0*Hr3(0,0,0) - 44.D0/3.D0*Hr3(0,0,1) - 44.D
     &    0/3.D0*Hr3(0,1,0) - 44.D0/3.D0*Hr3(1,0,0) )
      ppdiff = ppdiff + nf*cf**2 * ( 13.D0/9.D0*Hr1(0) - 13.D0/9.D0*
     &    Hr1(0)*x + 8.D0/9.D0*Hr2(0,0) - 28.D0/9.D0*Hr2(0,0)*x + 4.D0/
     &    3.D0*Hr2(0,1) - 4.D0/3.D0*Hr2(0,1)*x + 4.D0/3.D0*Hr2(1,0) - 4.
     &    D0/3.D0*Hr2(1,0)*x - Hr3(0,0,0) - Hr3(0,0,0)*x )
      ppdiff = ppdiff + nf*cf**2*pqq0p * (  - 11.D0/12.D0*Hr1(0) - 2.D0/
     &    3.D0*Hr1(0)*z2 + 11.D0/9.D0*Hr2(0,0) + 20.D0/9.D0*Hr2(0,1) + 
     &    20.D0/9.D0*Hr2(1,0) + 2.D0*Hr3(0,0,0) + 4.D0/3.D0*Hr3(0,0,1)
     &     + 4.D0/3.D0*Hr3(0,1,0) + 4.D0/3.D0*Hr3(1,0,0) )
*
       DTSP2NSP = 16.D0 * PPDIFF
*
       RETURN
       END
*
* ---------------------------------------------------------------------
*
*
* ..The third-order (NNLO) function P_NS-^T(2) - P_NS-^S(2) for the 
*   flavour differences of quark-antiquark differences and the total 
*   valence distribution.  
*
       FUNCTION DTSP2NSM (X, NF)
*
       IMPLICIT REAL*8 (A - Z)
       COMPLEX*16 HC1, HC2, HC3, HC4 
       INTEGER NF, N1, N2, NW, I1, I2, I3, N
       PARAMETER ( N1 = -1, N2 = 1, NW = 4 ) 
       DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
     ,           HC4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
     ,           HR4(N1:N2,N1:N2,N1:N2,N1:N2) 
       DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
     ,           HI4(N1:N2,N1:N2,N1:N2,N1:N2) 
       PARAMETER ( Z2 = 1.6449 34066 84822 64365 D0,
     ,             Z3 = 1.2020 56903 15959 42854 D0 )
*
* ...Colour factors
*
       CF  = 4./3.D0
       CA  = 3.D0
       CAM2CF = CA - 2.D0 * CF
*
* ...Some abbreviations
*
       PQQ0P = 2.D0/(1.D0-X) - 1.D0 - X
       PQQ0M = 2.D0/(1.D0+X) - 1.D0 + X
*
* ...The harmonic polylogs up to weight 4 by Gehrmann and Remiddi
*
       CALL HPLOG (X, NW, HC1,HC2,HC3,HC4, HR1,HR2,HR3,HR4,
     ,            HI1,HI2,HI3,HI4, N1, N2) 
*
* ...The splitting-function difference in terms of the harmonic polylogs
*
      pmdiff =
     &  + cf**2*cam2cf * (  - 178.D0/9.D0*Hr1(0) + 178.D0/9.D0*Hr1(0)
     &    *x - 3.D0*Hr1(0)*z3 + 3.D0*Hr1(0)*z3*x + 8.D0*Hr1(0)*z2 + 4.D0
     &    *Hr1(0)*z2*x - 140.D0/9.D0*Hr2(0,0) + 238.D0/9.D0*Hr2(0,0)*x
     &     - 4.D0*Hr2(0,0)*z2 - 52.D0/3.D0*Hr2(0,1) + 52.D0/3.D0*Hr2(0,
     &    1)*x - 52.D0/3.D0*Hr2(1,0) + 52.D0/3.D0*Hr2(1,0)*x + 8.D0*
     &    Hr3(-1,0,0) + 8.D0*Hr3(-1,0,0)*x + 4.D0*Hr3(0,-1,0) + 4.D0*
     &    Hr3(0,-1,0)*x - 13.D0/2.D0*Hr3(0,0,0) - 13.D0/2.D0*Hr3(0,0,0)
     &    *x - 8.D0*Hr3(0,0,1) - 8.D0*Hr3(0,0,1)*x - 4.D0*Hr3(0,1,0) - 
     &    4.D0*Hr3(0,1,0)*x - 4.D0*Hr4(0,-1,0,0) + 4.D0*Hr4(0,-1,0,0)*x
     &     - 4.D0*Hr4(0,0,-1,0) + 4.D0*Hr4(0,0,-1,0)*x + 8.D0*Hr4(0,0,0
     &    ,0) )
      pmdiff = pmdiff + cf**2*cam2cf*pqq0m * ( 7.D0*Hr1(0)*z3 - 3.D0/
     &    2.D0*Hr1(0)*z2 - 8.D0*Hr2(-1,0)*z2 - 8.D0*Hr2(0,-1)*z2 + 8.D0
     &    *Hr2(0,0)*z2 - 6.D0*Hr3(-1,0,0) - 3.D0*Hr3(0,-1,0) + 9.D0/2.D0
     &    *Hr3(0,0,0) - 16.D0*Hr4(-1,-1,0,0) - 8.D0*Hr4(-1,0,-1,0) + 18.
     &    D0*Hr4(-1,0,0,0) + 8.D0*Hr4(-1,0,0,1) + 4.D0*Hr4(-1,0,1,0) - 
     &    8.D0*Hr4(0,-1,-1,0) + 14.D0*Hr4(0,-1,0,0) + 4.D0*Hr4(0,-1,0,1
     &    ) + 8.D0*Hr4(0,0,-1,0) - 8.D0*Hr4(0,0,0,0) - 6.D0*Hr4(0,0,0,1
     &    ) - 2.D0*Hr4(0,0,1,0) )
      pmdiff = pmdiff + cf**2*cam2cf*pqq0p * ( 151.D0/24.D0*Hr1(0) + 
     &    Hr1(0)*z3 + 13.D0/6.D0*Hr1(0)*z2 - 169.D0/18.D0*Hr2(0,0) + 8.D
     &    0*Hr2(0,0)*z2 - 134.D0/9.D0*Hr2(0,1) + 4.D0*Hr2(0,1)*z2 - 134.
     &    D0/9.D0*Hr2(1,0) + 4.D0*Hr2(1,0)*z2 - 13.D0/2.D0*Hr3(0,0,0)
     &     - 22.D0/3.D0*Hr3(0,0,1) - 22.D0/3.D0*Hr3(0,1,0) - 22.D0/3.D0
     &    *Hr3(1,0,0) - 8.D0*Hr4(0,0,0,0) - 6.D0*Hr4(0,0,0,1) - 2.D0*
     &    Hr4(0,0,1,0) - 2.D0*Hr4(0,1,0,0) - 6.D0*Hr4(1,0,0,0) )
      pmdiff = pmdiff + cf**3 * (  - 325.D0/18.D0*Hr1(0) + 325.D0/18.D0
     &    *Hr1(0)*x + 3.D0*Hr1(0)*z2 - 5.D0*Hr1(0)*z2*x - 173.D0/18.D0*
     &    Hr2(0,0) + 691.D0/18.D0*Hr2(0,0)*x - 4.D0*Hr2(0,0)*z2 - 4.D0*
     &    Hr2(0,0)*z2*x - 50.D0/3.D0*Hr2(0,1) + 50.D0/3.D0*Hr2(0,1)*x
     &     - 50.D0/3.D0*Hr2(1,0) + 50.D0/3.D0*Hr2(1,0)*x + 25.D0/2.D0*
     &    Hr3(0,0,0) + 25.D0/2.D0*Hr3(0,0,0)*x + 2.D0*Hr3(0,0,1) + 2.D0
     &    *Hr3(0,0,1)*x + Hr3(0,1,0) + Hr3(0,1,0)*x )
      pmdiff = pmdiff + cf**3*pqq0p * ( 311.D0/24.D0*Hr1(0) + 4.D0/3.D0
     &    *Hr1(0)*z2 - 169.D0/9.D0*Hr2(0,0) + 8.D0*Hr2(0,0)*z2 - 268.D0/
     &    9.D0*Hr2(0,1) + 8.D0*Hr2(0,1)*z2 - 268.D0/9.D0*Hr2(1,0) + 8.D0
     &    *Hr2(1,0)*z2 - 22.D0*Hr3(0,0,0) - 44.D0/3.D0*Hr3(0,0,1) - 44.D
     &    0/3.D0*Hr3(0,1,0) - 44.D0/3.D0*Hr3(1,0,0) )
      pmdiff = pmdiff + nf*cf**2 * ( 13.D0/9.D0*Hr1(0) - 13.D0/9.D0*
     &    Hr1(0)*x + 8.D0/9.D0*Hr2(0,0) - 28.D0/9.D0*Hr2(0,0)*x + 4.D0/
     &    3.D0*Hr2(0,1) - 4.D0/3.D0*Hr2(0,1)*x + 4.D0/3.D0*Hr2(1,0) - 4.
     &    D0/3.D0*Hr2(1,0)*x - Hr3(0,0,0) - Hr3(0,0,0)*x )
      pmdiff = pmdiff + nf*cf**2*pqq0p * (  - 11.D0/12.D0*Hr1(0) - 2.D0/
     &    3.D0*Hr1(0)*z2 + 11.D0/9.D0*Hr2(0,0) + 20.D0/9.D0*Hr2(0,1) + 
     &    20.D0/9.D0*Hr2(1,0) + 2.D0*Hr3(0,0,0) + 4.D0/3.D0*Hr3(0,0,1)
     &     + 4.D0/3.D0*Hr3(0,1,0) + 4.D0/3.D0*Hr3(1,0,0) )

*
       DTSP2NSM = 16.D0 * PMDIFF
*
       RETURN
       END
*
* =================================================================av==
