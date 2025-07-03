
* ---------------------------------------------------------------------
*
*   Analytical results for the fourth-order coefficients Ag4, Bg4, 
*   Cg4 and Dg4 of the large-x expansion of the MSbar gluon-gluon 
*   splitting function P_gg(x) as defined and specified 
*   - including numerical values, see xPgg3a.f - in eqs. (13)-(21) of 
*
*       Additional moments and x-space approximations
*         of four-loop splitting functions in QCD
*       S. Moch, B. Ruijl, T. Ueda, J. Vermaseren and A. Vogt
*       DESY 23-150, Nikhef 2023-016, LTH 1354
*
*   Ag4 is identical to the four-loop cusp anomalous dimension, 
*   already determined by others in refs. [29,30] of the paper. 
*   The nf parts of Bg4 were given before in 
*
*       G. Das, S. Moch and A. Vogt
*       arXiv:2004.00563 = Phys. Lett. B807 (2020) 135546,
*
*   the remaining results are new.
*
*   Notations: z`i'=zeta_`i' are values of Riemann's zeta-function.
*              nf denotes the number of light flavours, for the 
*              group invariants see the paper. The perturbative 
*              expansion is written in terms of alpha_s/(4*pi). 
* 
* ---------------------------------------------------------------------

 L  Ag4 =
       + [d4AA/na] * (  - 128*z2 + 128/3*z3 - 384*z3^2 + 3520/3*z5 - 992*z6 )

       + nf*[d4RA/na] * ( 256*z2 - 256/3*z3 - 1280/3*z5 )

       + ca*nf^3 * (  - 32/81 + 64/27*z3 )

       + ca*cf*nf^2 * ( 2392/81 - 640/9*z3 + 32*z4 )

       + ca*cf^2*nf * ( 572/9 + 592/3*z3 - 320*z5 )

       + ca^2*nf^2 * ( 923/81 - 608/81*z2 + 2240/27*z3 - 112/3*z4 )

       + ca^2*cf*nf * (  - 34066/81 + 440/3*z2 - 352/5*z2^2 + 3712/9*z3 
          - 128*z3*z2 + 160*z5 )

       + ca^3*nf * (  - 24137/81 + 20320/81*z2 - 352/15*z2^2 - 23104/27*z3 
          + 448/3*z3*z2 + 2096/9*z5 )

       + ca^4 * ( 84278/81 - 88400/81*z2 + 20944/27*z3 - 352/3*z3*z2 - 16*z3^2
          + 1804*z4 - 3608/9*z5 - 2504/3*z6 );

  L Bg4 =
       + [d4AA/na] * (  - 800/9 + [B4(D4RA/nc)] + 1184/3*z2 - 1016/15*z2^2 
          + 5984/315*z2^3 - 784/3*z3 - 272*z3*z2 + 760/3*z5 )

       + nf*[d4RA/na] * ( 1952/9 + [B4(nf*D4RR/nc)] - 2368/3*z2 + 1312/3*z3  
          + 544*z3*z2 + 1016/3*z4 - 1520/3*z5 - 1496/9*z6 )

       + nf^2 * (  - 704/9*[d4RR/na] + 512/3*[d4RR/na]*z3 )

       + cf*nf^3 * ( 154/243 )

       + cf^2*nf^2 * ( 338/27 - 176/9*z3 )

       + cf^3*nf * ( 23 )

       + ca*nf^3 * ( 5/243 )

       + ca*cf*nf^2 * ( 3910/243 + 160/9*z3 )

       + ca*cf^2*nf * (  - 2723/27 + [B4(nf*cf^3)] - 162*z2 + 2948/9*z3 
          + 256/3*z3*z2 - 224*z3^2 + 204*z4 - 912*z5 + 6434/9*z6 )

       + ca^2*nf^2 * ( 1352/81 + 37/27*z2 + 289/27*z3 - 32/9*z3*z2 + 200/27*z4
          - 8/9*z5 )

       + ca^2*cf*nf * ( 23566/243 + [B4(nf*cf^2*ca)] + 4198/27*z2 + 8854/27*z3
          - 2744/9*z3*z2 + 1928/3*z3^2 - 27269/27*z4 + 6712/9*z5 - 2879/9*z6 )

       + ca^3*nf * (  - 8075/108 - 1/48*[B4(nf*D4RR/nc)] 
          - 1/2*[B4(nf*cf^2*ca)] - 1/4*[B4(nf*cf^3)] - 6155/54*z2 
          - 22714/27*z3 + 1874/9*z3*z2 - 1268/3*z3^2 + 7789/18*z4 + 919/9*z5 
          + 1777/54*z6 )

       + ca^4 * ( 50387/486 - 1/24*[B4(D4RA/nc)] + 2098/27*z2 + 1793/27*z2^2
          - 76516/945*z2^3 + 48088/27*z3 - 3902/9*z3*z2 + 336/5*z3*z2^2 
          + 682/3*z3^2 - 14617/9*z5 + 80*z5*z2 + 700*z7 );

  L Cg4 =
       + ca^2*nf^2 * ( 1216/81 )

       + ca^2*cf*nf * (  - 880/3 + 256*z3 )

       + ca^3*nf * (  - 41504/81 + 640/3*z2 - 896/3*z3 )

       + ca^4 * ( 177664/81 - 4288/3*z2 + 1728/5*z2^2 + 704/3*z3 );

  L Dg4 = 
       + ca^2*nf^2 * (  - 64/27 )

       + ca^2*cf*nf * (  - 8 )

       + ca^3*nf * ( 2048/27 - 16*z2 - 16/3*z2^2 - 160*z3 )

       + ca^4 * (  - 1984/27 + 16*z2 + 88/3*z2^2 + 1072*z3 - 160*z3*z2 
          - 320*z5 );

* -----------------------------------------------------------------------------
