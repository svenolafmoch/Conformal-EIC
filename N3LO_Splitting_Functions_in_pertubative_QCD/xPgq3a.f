*
* ..File: xPgq3a.f   
*
* ..The 4-loop MSbar quark-to-gluon splitting function P_{gq}^(3) 
*    for the evolution of the flavour-singlet parton distributions.
*    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
*
* ..These are approximations for fixed nf = 3, 4 and 5 based on the 
*    first five even moments together with small-x/large-x constraints.
*    The two sets providing the error estimate are called via IMOD = 1 
*    and IMOD = 2.  Any other value of IMOD invokes their average.
*
*    Reference: Additional moments and x-space approximations
*                 of four-loop splitting functions in QCD
*               S. Moch, B. Ruijl, T. Ueda, J. Vermaseren and A. Vogt
*               DESY 23-150, Nikhef 2023-016, LTH 1354   
*
* =====================================================================
*
*
       FUNCTION P3GQA (Y, NF, IMOD)
*
       IMPLICIT REAL*8 (A-Z)
       INTEGER IMOD, nf,nf2,nf3
*
       YM  = 1.D0/Y
       Y1  = 1.D0-Y
       DL  = DLOG(Y)
       DL1 = DLOG(1.D0-Y)
*
       nf2 = nf*nf
       nf3 = nf*nf2
*
* Known large-x coefficients
       x1L5cff = 1.3443073E+1 - 5.4869684E-1*nf
       x1L4cff = 3.7539831E+2 - 3.4494742E+1*nf + 8.7791495E-1*nf2
*    
* Small-x, Casimir scaled from P_gg (approx. for bfkl1)
       bfkl0 =   - 8.3086173E+3 / 2.25
       bfkl1 = ( - 1.0691199E+5 - nf* 9.9638304E+2 ) / 2.25
*
* The resulting part of the function
*
       P3gq01 =
     ,      + bfkl0* DL**3*YM
     ,      + bfkl1* DL**2*YM
     ,      + x1L4cff* DL1**4
     ,      + x1L5cff* DL1**5
*
* The selected approximations for nf = 3, 4, 5
       IF ( NF .EQ. 3 ) THEN
         P3gqApp1 = P3GQ01
     ,      + 3.4*bfkl1* DL*YM
     ,      - 161562.* Y1*YM
     ,      + 36469.
     ,      + 72317.* DL
     ,      - 3977.3* DL1**2
     ,      +  484.4* DL1**3
         P3gqApp2  = P3GQ01
     ,      + 5.4*bfkl1* DL*YM
     ,      - 546482.* Y1*YM
     ,      - 39464.
     ,      - 401000.* DL
     ,      + 13270.* DL1**2
     ,      + 3289.* DL1**3
       ELSE IF ( NF .EQ. 4 ) THEN
         P3gqApp1 = P3GQ01
     ,      + 3.4*bfkl1* DL*YM
     ,      - 158805.* Y1*YM
     ,      + 35098.
     ,      + 87258.* DL
     ,      - 4834.1* DL1**2
     ,      +  176.6* DL1**3
         P3gqApp2  = P3GQ01
     ,      + 5.4*bfkl1* DL*YM
     ,      - 547215.* Y1*YM
     ,      - 41523.
     ,      - 390350.* DL
     ,      + 12571.* DL1**2
     ,      + 3007.* DL1**3
       ELSE IF ( NF .EQ. 5 ) THEN
         P3gqApp1 = P3GQ01
     ,      + 3.4*bfkl1* DL*YM
     ,      - 154336.* Y1*YM
     ,      + 33889.  
     ,      + 103440.* DL
     ,      - 5745.8 * DL1**2
     ,      -  128.6* DL1**3
         P3gqApp2  = P3GQ01
     ,      + 5.4*bfkl1* DL*YM
     ,      - 546236.* Y1*YM
     ,      - 43421.
     ,      - 378460.* DL
     ,      + 11816.* DL1**2
     ,      + 2727.3* DL1**3
       ELSE IF ( NF .EQ. 5 ) THEN
         P3gqApp1 = P3GQ01
         P3gqApp2 = P3GQ01
       ELSE
         WRITE(6,*) '  Error in function P3GQA: choice of nf   '
         CALL ABORT
       END IF
*
* We return (for now) one of the two error-band representatives
* or the present best estimate, their average
       IF ( IMOD .EQ. 1 ) THEN
         P3GQA = P3gqApp1
       ELSE IF ( IMOD .EQ. 2 ) THEN
         P3GQA = P3gqApp2
       ELSE
         P3GQA = 0.5* ( P3gqApp1 + P3gqApp2 )
       END IF
*
       RETURN
       END
*
* =================================================================av==
