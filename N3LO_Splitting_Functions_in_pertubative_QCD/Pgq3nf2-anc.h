*----------------------------------------------------------------------*
*
*  The nf^2 contribution to the 4-loop splitting function P_gq^(3)
*  for the N^3LO evolution of unpolarized parton distributions 
*
*  H(R(m1,...,mk),x) denote harmonic polylogarithms, with the indices
*  abbreviated as in eq. (7) og the paper, H(R(2),x) = H(R(0,1),x) etc;
*  z2, ..., z5 are values of the Riemann zeta-function.
*
*  Besides the all-x result, also the expansions in the small-x limit 
*  (x^{-1,0} terms) and the large-x limit ([1-x]^0 terms) are given.
*
*  Reference: G. Falcioni, F. Herzog, S. Moch, J. Vermaseren and A. Vogt
*               The double fermionic contribution to the
*               four-loop quark-to-gluon splitting function
*             DESY 23-146, LTH 1353
*
*----------------------------------------------------------------------*

 L  Pgq3nf2 = 16/27* cf^2*nf^2 * (

       + H(R(3,0,0),x) * (
          + 48
          - 24*x
          )
       + H(R(3,1,0),x) * (
          + 24
          - 12*x
          )
       + H(R(3,1,1),x) * (
          + 24
          - 12*x
          )
       + H(R(3,2),x) * (
          + 24
          - 12*x
          )
       + H(R(4,0),x) * (
          + 120
          - 60*x
          )
       + H(R(5),x) * (
          + 312
          - 156*x
          )
       + H(R(-3,0),x) * (
          - 96
          )
       + H(R(-2,-1,0),x) * (
          + 288
          - 96*x^-1
          - 144*x
          )
       + H(R(-2,0,0),x) * (
          - 192
          + 64*x^-1
          + 96*x
          )
       + H(R(-2,2),x) * (
          - 288
          + 96*x^-1
          + 144*x
          )
       + H(R(0,0,0,0),x) * (
          - 766
          + 485*x
          )
       + H(R(1,0,0,0),x) * (
          + 16
          - 16*x^-1
          - 8*x
          )
       + H(R(1,1,0,0),x) * (
          - 32
          + 32*x^-1
          + 16*x
          )
       + H(R(1,1,1,0),x) * (
          - 72
          + 72*x^-1
          + 36*x
          )
       + H(R(1,1,1,1),x) * (
          + 32
          - 32*x^-1
          - 16*x
          )
       + H(R(1,1,2),x) * (
          - 56
          + 56*x^-1
          + 28*x
          )
       + H(R(1,2,0),x) * (
          - 16
          + 16*x^-1
          + 8*x
          )
       + H(R(1,2,1),x) * (
          - 32
          + 32*x^-1
          + 16*x
          )
       + H(R(1,3),x) * (
          - 16
          + 16*x^-1
          + 8*x
          )
       + H(R(2,0,0),x) * (
          - 100
          - 64*x^-1
          - 10*x
          )
       + H(R(2,1,0),x) * (
          - 20
          - 16*x^-1
          + 28*x
          )
       + H(R(2,1,1),x) * (
          - 22
          - 16*x^-1
          + 29*x
          )
       + H(R(2,2),x) * (
          - 20
          - 16*x^-1
          + 28*x
          )
       + H(R(3,0),x) * (
          - 80
          + 130*x
          )
       + H(R(3,1),x) * (
          - 102
          + 81*x
          )
       + H(R(4),x) * (
          + 14
          - 31*x
          )
       + H(R(-2,0),x) * (
          + 480
          - 424/3*x^-1
          - 312*x
          + 64/3*x^2
          )
       + H(R(-1,-1,0),x) * (
          - 248*x^-1
          + 216*x
          - 32*x^2
          )
       + H(R(-1,0,0),x) * (
          + 496/3*x^-1
          - 144*x
          + 64/3*x^2
          )
       + H(R(-1,2),x) * (
          + 248*x^-1
          - 216*x
          + 32*x^2
          )
       + H(R(0,0,0),x) * (
          + 563/3
          - 91/3*x
          + 160/3*x^2
          - 312*z2
          + 156*z2*x
          )
       + H(R(1,0,0),x) * (
          + 764/3
          - 244*x^-1
          - 259/3*x
          + 64/3*x^2
          )
       + H(R(1,1,0),x) * (
          + 542/3
          - 152*x^-1
          - 361/3*x
          + 16/3*x^2
          )
       + H(R(1,1,1),x) * (
          - 268/3
          + 116*x^-1
          + 409/6*x
          + 16/3*x^2
          )
       + H(R(1,2),x) * (
          + 418/3
          - 332/3*x^-1
          - 275/3*x
          + 16/3*x^2
          )
       + H(R(2,0),x) * (
          + 1552/3
          + 368/3*x^-1
          - 149/3*x
          + 112/3*x^2
          )
       + H(R(2,1),x) * (
          + 161/3
          + 36*x^-1
          - 244/3*x
          + 32/3*x^2
          )
       + H(R(3),x) * (
          + 1335
          + 390*x
          + 160/3*x^2
          - 24*z2
          + 12*z2*x
          )
       + H(R(-2),x) * (
          - 144*z2*x^-1
          + 432*z2
          - 216*z2*x
          )
       + H(R(-1,0),x) * (
          - 160
          - 2408/9*x^-1
          + 56*x
          - 464/9*x^2
          )
       + H(R(0,0),x) * (
          - 8110/3
          - 4649/6*x
          - 2080/9*x^2
          - 14*z2
          + 31*z2*x
          - 72*z3
          + 36*z3*x
          )
       + H(R(1,0),x) * (
          - 850/3
          + 3776/9*x^-1
          + 41/6*x
          - 896/9*x^2
          - 16*z2*x^-1
          + 16*z2
          - 8*z2*x
          )
       + H(R(1,1),x) * (
          - 3461/18
          + 491/9*x^-1
          + 2027/36*x
          - 376/9*x^2
          - 56*z2*x^-1
          + 56*z2
          - 28*z2*x
          )
       + H(R(2),x) * (
          - 6727/9
          - 3056/9*x^-1
          - 2891/18*x
          - 1384/9*x^2
          + 64*z2*x^-1
          + 164*z2
          + 44*z2*x
          )
       + H(R(-1),x) * (
          - 372*z2*x^-1
          + 324*z2*x
          - 48*z2*x^2
          )
       + H(R(0),x) * (
          + 1344
          + 10528/27*x^-1
          + 389/36*x
          + 6494/27*x^2
          - 424/3*z2*x^-1
          - 1335*z2
          - 702*z2*x
          - 160/3*z2*x^2
          + 96*z3*x^-1
          - 48*z3
          + 342*z3*x
          - 48*z4
          + 24*z4*x
          )
       + H(R(1),x) * (
          + 72371/36
          - 93103/54*x^-1
          - 30551/72*x
          + 2678/27*x^2
          + 704/3*z2*x^-1
          - 418/3*z2
          - 49/3*z2*x
          - 64/3*z2*x^2
          + 32*z3*x^-1
          - 32*z3
          + 16*z3*x
          )
       - 363779/216
          + 368893/324*x^-1
          + 326947/432*x
          - 7090/81*x^2
          + 72*z2*x^-1
          + 6727/9*z2
          + 3899/18*z2*x
          + 1384/9*z2*x^2
          - 116/3*z3*x^-1
          + 674*z3
          - 1288*z3*x
          + 80/3*z3*x^2
          - 24*z3*z2
          + 12*z3*z2*x
          + 148*z4*x^-1
          - 1229/2*z4
          + 1295/4*z4*x
          - 24*z5
          + 12*z5*x
        )

     + 16/27* cf*ca*nf^2 * (

       + H(R(-3,0),x) * (
          + 80
          )
       + H(R(-2,-1,0),x) * (
          - 192
          + 48*x^-1
          + 72*x
          )
       + H(R(-2,0,0),x) * (
          + 112
          - 32*x^-1
          - 48*x
          )
       + H(R(-2,2),x) * (
          + 192
          - 48*x^-1
          - 72*x
          )
       + H(R(-1,-2,0),x) * (
          - 24
          - 24*x^-1
          - 12*x
          )
       + H(R(-1,-1,-1,0),x) * (
          + 116
          + 116*x^-1
          + 58*x
          )
       + H(R(-1,-1,0,0),x) * (
          - 48
          - 48*x^-1
          - 24*x
          )
       + H(R(-1,-1,2),x) * (
          - 92
          - 92*x^-1
          - 46*x
          )
       + H(R(-1,0,0,0),x) * (
          - 4
          - 4*x^-1
          - 2*x
          )
       + H(R(-1,2,0),x) * (
          + 4
          + 4*x^-1
          + 2*x
          )
       + H(R(-1,2,1),x) * (
          + 40
          + 40*x^-1
          + 20*x
          )
       + H(R(-1,3),x) * (
          + 46
          + 46*x^-1
          + 23*x
          )
       + H(R(0,0,0,0),x) * (
          + 12
          - 6*x
          )
       + H(R(1,-2,0),x) * (
          + 40
          - 40*x^-1
          - 20*x
          )
       + H(R(1,0,0,0),x) * (
          + 4
          - 4*x^-1
          - 2*x
          )
       + H(R(1,1,0,0),x) * (
          + 28
          - 28*x^-1
          - 14*x
          )
       + H(R(1,1,1,0),x) * (
          - 54
          + 54*x^-1
          + 27*x
          )
       + H(R(1,1,1,1),x) * (
          - 32
          + 32*x^-1
          + 16*x
          )
       + H(R(1,1,2),x) * (
          - 50
          + 50*x^-1
          + 25*x
          )
       + H(R(1,2,0),x) * (
          - 16
          + 16*x^-1
          + 8*x
          )
       + H(R(1,2,1),x) * (
          - 52
          + 52*x^-1
          + 26*x
          )
       + H(R(1,3),x) * (
          + 50
          - 50*x^-1
          - 25*x
          )
       + H(R(2,0,0),x) * (
          + 112
          + 32*x^-1
          + 44*x
          )
       + H(R(2,1,0),x) * (
          - 12
          + 24*x^-1
          + 22*x
          )
       + H(R(2,1,1),x) * (
          - 140
          + 52*x^-1
          + 4*x
          )
       + H(R(2,2),x) * (
          - 28
          + 16*x^-1
          - 2*x
          )
       + H(R(3,0),x) * (
          + 24
          + 2*x
          )
       + H(R(3,1),x) * (
          + 36
          + 22*x
          )
       + H(R(4),x) * (
          + 36
          - 16*x
          )
       + H(R(-2,0),x) * (
          - 956/3
          + 212/3*x^-1
          + 152*x
          - 32/3*x^2
          )
       + H(R(-1,-1,0),x) * (
          + 460/3
          + 958/3*x^-1
          - 34/3*x
          + 16*x^2
          )
       + H(R(-1,0,0),x) * (
          - 260/3
          - 550/3*x^-1
          + 50/3*x
          - 32/3*x^2
          )
       + H(R(-1,2),x) * (
          - 140/3
          - 638/3*x^-1
          + 242/3*x
          - 16*x^2
          )
       + H(R(0,0,0),x) * (
          - 332/3
          - 19/12*x
          + 8/3*x^2
          )
       + H(R(1,0,0),x) * (
          - 184/3
          + 434/3*x^-1
          - 106/3*x
          - 40/3*x^2
          )
       + H(R(1,1,0),x) * (
          + 283/3
          - 83*x^-1
          - 248/3*x
          - 4/3*x^2
          )
       + H(R(1,1,1),x) * (
          + 700/3
          - 259*x^-1
          - 299/3*x
          + 44/3*x^2
          )
       + H(R(1,2),x) * (
          + 277/3
          - 311/3*x^-1
          - 188/3*x
          + 4/3*x^2
          )
       + H(R(2,0),x) * (
          + 67/3
          - 212/3*x^-1
          - 265/6*x
          - 4/3*x^2
          )
       + H(R(2,1),x) * (
          + 525/2
          - 76*x^-1
          - 7/4*x
          - 12*x^2
          )
       + H(R(3),x) * (
          - 1499/6
          - 343/2*x
          + 4/3*x^2
          )
       + H(R(-2),x) * (
          + 72*z2*x^-1
          - 288*z2
          + 108*z2*x
          )
       + H(R(-1,-1),x) * (
          + 150*z2*x^-1
          + 150*z2
          + 75*z2*x
          )
       + H(R(-1,0),x) * (
          + 4/3*x^-2
          + 226/3
          + 1498/9*x^-1
          + 7/3*x
          + 232/9*x^2
          - 38*z2*x^-1
          - 38*z2
          - 19*z2*x
          )
       + H(R(0,0),x) * (
          + 8735/18
          + 3091/36*x
          + 440/9*x^2
          - 36*z2
          + 16*z2*x
          )
       + H(R(1,0),x) * (
          + 1363/18
          - 305/6*x^-1
          + 301/18*x
          + 52/9*x^2
          + 42*z2*x^-1
          - 42*z2
          + 21*z2*x
          )
       + H(R(1,1),x) * (
          - 8747/36
          + 12383/36*x^-1
          + 1498/9*x
          - 10*x^2
          + 8*z2*x^-1
          - 8*z2
          + 4*z2*x
          )
       + H(R(2),x) * (
          + 16495/36
          + 1067/6*x^-1
          + 610/9*x
          + 1202/9*x^2
          - 40*z2*x^-1
          - 68*z2
          - 34*z2*x
          )
       + H(R(-1),x) * (
          + 1117/3*z2*x^-1
          + 370/3*z2
          - 259/3*z2*x
          + 24*z2*x^2
          - 150*z3*x^-1
          - 150*z3
          - 75*z3*x
          )
       + H(R(0),x) * (
          - 10607/12
          - 4916/27*x^-1
          - 12139/72*x
          - 9944/27*x^2
          + 212/3*z2*x^-1
          + 1499/6*z2
          + 647/2*z2*x
          - 4/3*z2*x^2
          - 72*z3*x^-1
          - 72*z3
          - 146*z3*x
          )
       + H(R(1),x) * (
          - 4495/6
          + 132025/216*x^-1
          + 3983/36*x
          - 6500/27*x^2
          - 56*z2*x^-1
          - 47/3*z2
          + 205/3*z2*x
          + 20/3*z2*x^2
          - 66*z3*x^-1
          + 66*z3
          - 33*z3*x
          )
       + 30829/108
          - 191575/648*x^-1
          - 91561/432*x
          + 27377/81*x^2
          - 205/18*z2*x^-1
          - 16495/36*z2
          - 589/9*z2*x
          - 1202/9*z2*x^2
          + 104*z3*x^-1
          - 752*z3
          + 8587/12*z3*x
          + 8*z3*x^2
          - 114*z4*x^-1
          + 331*z4
          - 59/2*z4*x
     );


 L  Pgq3nf2xto0 =

       + 1/x*ln_(x) * (
          - 78080/729*nf^2*cf*ca
          + 3392/81*nf^2*cf*ca*z2
          - 128/3*nf^2*cf*ca*z3
          + 168448/729*nf^2*cf^2
          - 6784/81*nf^2*cf^2*z2
          + 512/9*nf^2*cf^2*z3
          )
       + 1/x * (
          - 384878/2187*nf^2*cf*ca
          - 1640/243*nf^2*cf*ca*z2
          + 1664/27*nf^2*cf*ca*z3
          - 608/9*nf^2*cf*ca*z4
          + 1475572/2187*nf^2*cf^2
          + 128/3*nf^2*cf^2*z2
          - 1856/81*nf^2*cf^2*z3
          + 2368/27*nf^2*cf^2*z4
          )
       + ln_(x)^4 * (
          + 8/27*nf^2*cf*ca
          - 1532/81*nf^2*cf^2
          )
       + ln_(x)^3 * (
          - 2848/243*nf^2*cf*ca
          + 4120/243*nf^2*cf^2
          - 832/27*nf^2*cf^2*z2
          )
       + ln_(x)^2 * (
          + 32732/243*nf^2*cf*ca
          - 32/3*nf^2*cf*ca*z2
          - 7376/9*nf^2*cf^2
          - 112/27*nf^2*cf^2*z2
          - 64/3*nf^2*cf^2*z3
          )
       + ln_(x) * (
          - 106316/243*nf^2*cf*ca
          + 12184/81*nf^2*cf*ca*z2
          - 128/3*nf^2*cf*ca*z3
          + 24640/27*nf^2*cf^2
          - 21616/27*nf^2*cf^2*z2
          - 256/9*nf^2*cf^2*z3
          - 256/9*nf^2*cf^2*z4
          )
       + 133714/243*nf^2*cf*ca
          - 16396/243*nf^2*cf*ca*z2
          - 15488/27*nf^2*cf*ca*z3
          + 5296/27*nf^2*cf*ca*z4
          - 565226/243*nf^2*cf^2
          + 78640/243*nf^2*cf^2*z2
          + 11296/27*nf^2*cf^2*z3
          - 128/9*nf^2*cf^2*z3*z2
          - 9832/27*nf^2*cf^2*z4
          - 128/9*nf^2*cf^2*z5
         ;

 L  Pgq3nf2xto1 =

      +ln_(1-x)^4*(
         +32/81*cf*ca*nf^2
         -32/81*cf^2*nf^2
         )
      +ln_(1-x)^3*(
         +2656/243*cf*ca*nf^2
         -2404/243*cf^2*nf^2
         )
      +ln_(1-x)^2*(
         +18536/243*cf*ca*nf^2
         +16/27*z2*cf*ca*nf^2
         -8870/243*cf^2*nf^2
         -32/3*z2*cf^2*nf^2
         )
      +ln_(1-x)*(
         +12866/81*cf*ca*nf^2
         -160/81*z2*cf*ca*nf^2
         +880/27*z3*cf*ca*nf^2
         +1870/81*cf^2*nf^2
         -4144/81*z2*cf^2*nf^2
         -320/27*z3*cf^2*nf^2
         )
      +5561/81*cf*ca*nf^2
         -2048/243*z2*cf*ca*nf^2
         +12800/81*z3*cf*ca*nf^2
         -3664/135*z2^2*cf*ca*nf^2
         +1979/27*cf^2*nf^2
         -232/9*z2*cf^2*nf^2
         -9472/81*z3*cf^2*nf^2
         +640/27*z2^2*cf^2*nf^2
         ;

*----------------------------------------------------------------------*
