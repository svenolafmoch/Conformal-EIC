* ----------------------------------------------------------------------
*
*  The nf^2 contribution to the N^3LO anomalous dimension gamma_gq^(3)
*  for the evolution of unpolarized parton distributions at all even N.
* 
*  S(R(m1,...,mk),N) denote harmonic sums, den(x) = 1/x denominators,
*  z3 and z4 are values ot Riemann's zeta function.
*
*  Reference: G. Falcioni, F. Herzog, S. Moch, J. Vermaseren and A. Vogt
*               The double fermionic contribution to the 
*               four-loop quark-to-gluon splitting function
*             DESY 23-146, LTH 1353
*
* ---------------------------------------------------------------------- 

 L  ggq3nf2 = 16/27* cf^2*nf^2* (

         S(R(1,1,1,1),N) * (
          + 32*den( - 1 + N)
          + 16*den(1 + N)
          - 32*den(N)
          )
       + S(R(1,1,2),N) * (
          + 72*den( - 1 + N)
          + 36*den(1 + N)
          - 72*den(N)
          )
       + S(R(1,2,1),N) * (
          + 56*den( - 1 + N)
          + 28*den(1 + N)
          - 56*den(N)
          )
       + S(R(1,3),N) * (
          - 32*den( - 1 + N)
          - 16*den(1 + N)
          + 32*den(N)
          )
       + S(R(2,1,1),N) * (
          + 32*den( - 1 + N)
          + 16*den(1 + N)
          - 32*den(N)
          )
       + S(R(2,2),N) * (
          - 16*den( - 1 + N)
          - 8*den(1 + N)
          + 16*den(N)
          )
       + S(R(3,1),N) * (
          - 16*den( - 1 + N)
          - 8*den(1 + N)
          + 16*den(N)
          )
       + S(R(4),N) * (
          - 16*den( - 1 + N)
          - 8*den(1 + N)
          + 16*den(N)
          )
       + S(R(-3),N) * (
          - 496/3*den( - 1 + N)
          + 64*den( - 1 + N)^2
          + 144*den(1 + N)
          + 96*den(1 + N)^2
          + 64/3*den(2 + N)
          + 192*den(N)^2
          )
       + S(R(-2,1),N) * (
          + 248*den( - 1 + N)
          - 96*den( - 1 + N)^2
          - 216*den(1 + N)
          - 144*den(1 + N)^2
          - 32*den(2 + N)
          - 288*den(N)^2
          )
       + S(R(1,-2),N) * (
          + 248*den( - 1 + N)
          - 96*den( - 1 + N)^2
          - 216*den(1 + N)
          - 144*den(1 + N)^2
          - 32*den(2 + N)
          - 288*den(N)^2
          )
       + S(R(1,1,1),N) * (
          - 148*den( - 1 + N)
          - 16*den( - 1 + N)^2
          - 409/6*den(1 + N)
          + 45*den(1 + N)^2
          + 12*den(1 + N)^3
          - 16/3*den(2 + N)
          + 364/3*den(N)
          - 22*den(N)^2
          - 24*den(N)^3
          )
       + S(R(1,2),N) * (
          - 224*den( - 1 + N)
          + 16*den( - 1 + N)^2
          - 361/3*den(1 + N)
          + 8*den(1 + N)^2
          - 12*den(1 + N)^3
          + 16/3*den(2 + N)
          + 758/3*den(N)
          + 20*den(N)^2
          + 24*den(N)^3
          )
       + S(R(2,1),N) * (
          - 500/3*den( - 1 + N)
          + 16*den( - 1 + N)^2
          - 275/3*den(1 + N)
          - 12*den(1 + N)^3
          + 16/3*den(2 + N)
          + 586/3*den(N)
          + 20*den(N)^2
          + 24*den(N)^3
          )
       + S(R(3),N) * (
          + 276*den( - 1 + N)
          - 64*den( - 1 + N)^2
          + 259/3*den(1 + N)
          - 26*den(1 + N)^2
          + 24*den(1 + N)^3
          - 64/3*den(2 + N)
          - 860/3*den(N)
          - 100*den(N)^2
          - 48*den(N)^3
          )
       + S(R(-2),N) * (
          - 5504/9*den( - 1 + N)
          + 712/3*den( - 1 + N)^2
          + 24*den(1 + N)
          + 96*den(1 + N)^2
          - 144*den(1 + N)^3
          + 752/9*den(2 + N)
          - 32/3*den(2 + N)^2
          + 504*den(N)
          + 480*den(N)^2
          + 96*den(N)^3
          )
       + S(R(1,1),N) * (
          + 121/9*den( - 1 + N)
          + 52*den( - 1 + N)^2
          - 2219/36*den(1 + N)
          - 299/2*den(1 + N)^2
          - 20*den(1 + N)^3
          + 12*den(1 + N)^4
          + 424/9*den(2 + N)
          + 16/3*den(2 + N)^2
          + 2237/18*den(N)
          + 257/3*den(N)^2
          + 102*den(N)^3
          )
       + S(R(2),N) * (
          + 5432/9*den( - 1 + N)
          - 416/3*den( - 1 + N)^2
          + 73/6*den(1 + N)
          - 212/3*den(1 + N)^2
          + 130*den(1 + N)^3
          + 48*den(1 + N)^4
          - 944/9*den(2 + N)
          - 32*den(2 + N)^2
          - 1402/3*den(N)
          - 1600/3*den(N)^2
          - 80*den(N)^3
          - 120*den(N)^4
          )
       + S(R(1),N) * (
          + 82801/54*den( - 1 + N)
          - 2660/9*den( - 1 + N)^2
          + 31639/72*den(1 + N)
          - 2219/12*den(1 + N)^2
          - 2491/6*den(1 + N)^3
          + 85*den(1 + N)^4
          + 156*den(1 + N)^5
          - 3086/27*den(2 + N)
          - 352/3*den(2 + N)^2
          - 224/3*den(2 + N)^3
          - 65503/36*den(N)
          - 4915/9*den(N)^2
          - 1351*den(N)^3
          + 14*den(N)^4
          - 312*den(N)^5
          )
       - 1094983/324*den( - 1 + N)
          + 19192/27*den( - 1 + N)^2
          - 125953/144*den(1 + N)
          + 2265/8*den(1 + N)^2
          + 6425/12*den(1 + N)^3
          - 4019/6*den(1 + N)^4
          - 256*den(1 + N)^5
          + 228*den(1 + N)^6
          + 16636/81*den(2 + N)
          + 312*den(2 + N)^2
          + 1312/9*den(2 + N)^3
          - 64*den(2 + N)^4
          + 282613/72*den(N)
          + 1540*den(N)^2
          + 2766*den(N)^3
          + 515/3*den(N)^4
          + 766*den(N)^5

     + 4*z3* (
         S(R(1),N) * (
          - 12*den( - 1 + N)
          - 6*den(1 + N)
          + 12*den(N)
          )
       + 52*den( - 1 + N)
          + 12*den( - 1 + N)^2
          + 76*den(1 + N)
          + 9*den(1 + N)^2
          + 4*den(2 + N)
          - 89*den(N)
          + 42*den(N)^2
          )

     + z4* (
          - 144*den( - 1 + N)
          + 144*den(N)
          - 72*den(1 + N)
          )
     )

     + 16/27* cf*ca*nf^2* (

       + S(R(-4),N) * (
          - 4*den( - 1 + N)
          - 2*den(1 + N)
          + 4*den(N)
          )
       + S(R(-3,1),N) * (
          - 46*den( - 1 + N)
          - 23*den(1 + N)
          + 46*den(N)
          )
       + S(R(-2,-2),N) * (
          - 40*den( - 1 + N)
          - 20*den(1 + N)
          + 40*den(N)
          )
       + S(R(-2,1,1),N) * (
          + 40*den( - 1 + N)
          + 20*den(1 + N)
          - 40*den(N)
          )
       + S(R(-2,2),N) * (
          - 4*den( - 1 + N)
          - 2*den(1 + N)
          + 4*den(N)
          )
       + S(R(1,-3),N) * (
          - 48*den( - 1 + N)
          - 24*den(1 + N)
          + 48*den(N)
          )
       + S(R(1,-2,1),N) * (
          + 92*den( - 1 + N)
          + 46*den(1 + N)
          - 92*den(N)
          )
       + S(R(1,1,-2),N) * (
          + 116*den( - 1 + N)
          + 58*den(1 + N)
          - 116*den(N)
          )
       + S(R(1,1,1,1),N) * (
          - 32*den( - 1 + N)
          - 16*den(1 + N)
          + 32*den(N)
          )
       + S(R(1,1,2),N) * (
          + 54*den( - 1 + N)
          + 27*den(1 + N)
          - 54*den(N)
          )
       + S(R(1,2,1),N) * (
          + 50*den( - 1 + N)
          + 25*den(1 + N)
          - 50*den(N)
          )
       + S(R(1,3),N) * (
          + 28*den( - 1 + N)
          + 14*den(1 + N)
          - 28*den(N)
          )
       + S(R(2,-2),N) * (
          - 24*den( - 1 + N)
          - 12*den(1 + N)
          + 24*den(N)
          )
       + S(R(2,1,1),N) * (
          + 52*den( - 1 + N)
          + 26*den(1 + N)
          - 52*den(N)
          )
       + S(R(2,2),N) * (
          - 16*den( - 1 + N)
          - 8*den(1 + N)
          + 16*den(N)
          )
       + S(R(3,1),N) * (
          + 50*den( - 1 + N)
          + 25*den(1 + N)
          - 50*den(N)
          )
       + S(R(4),N) * (
          - 4*den( - 1 + N)
          - 2*den(1 + N)
          + 4*den(N)
          )
       + S(R(-3),N) * (
          + 694/3*den( - 1 + N)
          - 32*den( - 1 + N)^2
          - 50/3*den(1 + N)
          - 72*den(1 + N)^2
          - 32/3*den(2 + N)
          - 404/3*den(N)
          - 112*den(N)^2
          )
       + S(R(-2,1),N) * (
          - 914/3*den( - 1 + N)
          + 48*den( - 1 + N)^2
          + 242/3*den(1 + N)
          + 118*den(1 + N)^2
          + 16*den(2 + N)
          + 416/3*den(N)
          + 192*den(N)^2
          )
       + S(R(1,-2),N) * (
          - 1306/3*den( - 1 + N)
          + 48*den( - 1 + N)^2
          + 34/3*den(1 + N)
          + 130*den(1 + N)^2
          + 16*den(2 + N)
          + 808/3*den(N)
          + 192*den(N)^2
          )
       + S(R(1,1,1),N) * (
          + 291*den( - 1 + N)
          + 52*den( - 1 + N)^2
          + 299/3*den(1 + N)
          - 12*den(1 + N)^2
          - 44/3*den(2 + N)
          - 796/3*den(N)
          - 140*den(N)^2
          )
       + S(R(1,2),N) * (
          - 137*den( - 1 + N)
          - 24*den( - 1 + N)^2
          - 248/3*den(1 + N)
          + 5*den(1 + N)^2
          - 4/3*den(2 + N)
          + 445/3*den(N)
          + 12*den(N)^2
          )
       + S(R(2,1),N) * (
          - 461/3*den( - 1 + N)
          - 16*den( - 1 + N)^2
          - 188/3*den(1 + N)
          + 27*den(1 + N)^2
          + 4/3*den(2 + N)
          + 427/3*den(N)
          + 28*den(N)^2
          )
       + S(R(3),N) * (
          - 518/3*den( - 1 + N)
          + 32*den( - 1 + N)^2
          + 106/3*den(1 + N)
          + 58*den(1 + N)^2
          + 40/3*den(2 + N)
          + 268/3*den(N)
          + 112*den(N)^2
          )
       + S(R(-2),N) * (
          + 5380/9*den( - 1 + N)
          - 356/3*den( - 1 + N)^2
          + 55/3*den(1 + N)
          - 422/3*den(1 + N)^2
          + 138*den(1 + N)^3
          - 376/9*den(2 + N)
          + 16/3*den(2 + N)^2
          - 1520/3*den(N)
          - 1148/3*den(N)^2
          - 80*den(N)^3
          )
       + S(R(1,1),N) * (
          - 23147/36*den( - 1 + N)
          - 128*den( - 1 + N)^2
          - 1630/9*den(1 + N)
          + 1175/12*den(1 + N)^2
          - 28*den(1 + N)^3
          + 74/3*den(2 + N)
          - 80/3*den(2 + N)^2
          + 19511/36*den(N)
          + 709/2*den(N)^2
          - 36*den(N)^3
          )
       + S(R(2),N) * (
          + 169/6*den( - 1 + N)
          + 284/3*den( - 1 + N)^2
          + 277/18*den(1 + N)
          - 77/2*den(1 + N)^2
          + den(1 + N)^3
          + 64/9*den(2 + N)
          - 59/18*den(N)
          - 127/3*den(N)^2
          + 24*den(N)^3
          )
       + S(R(1),N) * (
          + 7145/216*den( - 1 + N)
          + 1331/6*den( - 1 + N)^2
          - 3527/36*den(1 + N)
          - 128*den(1 + N)^2
          + 1513/12*den(1 + N)^3
          - 87*den(1 + N)^4
          + 6158/27*den(2 + N)
          + 1532/9*den(2 + N)^2
          - 32/3*den(2 + N)^3
          + 3775/36*den(N)
          + 4099/36*den(N)^2
          + 1523/6*den(N)^3
          + 36*den(N)^4
          )
       + 304049/324*den( - 1 + N)
          - 19435/54*den( - 1 + N)^2
          + 28097/48*den(1 + N)
          - 17785/72*den(1 + N)^2
          - 6941/36*den(1 + N)^3
          + 836/3*den(1 + N)^4
          - 88*den(1 + N)^5
          - 57623/81*den(2 + N)
          - 2632/9*den(2 + N)^2
          + 256/3*den(2 + N)^3
          - 66833/72*den(N)
          - 26555/36*den(N)^2
          - 8183/18*den(N)^3
          - 356/3*den(N)^4
          - 12*den(N)^5

       + theta_(N-4)* S(R(-2), - 2 + N) * (
          - 4/3*den( - 2 + N)
          )

     + 6*z3* (

       + 1/3* delta_(N-2)
       + S(R(1),N) * (
          + 12*den( - 1 + N)
          - 12*den(N)
          + 6*den(1 + N)
          )
       + 63*den(N)
          - 54*den(1+N)
          - 103/3*den(-1+N)
          - 8/3*den(2+N)
          - 28*den(N)^2
          - 12*den(1+N)^2
          - 8*den(-1+N)^2
          )

     + z4* (
          + 144*den( - 1 + N)
          - 144*den(N)
          + 72*den(1 + N)
          )
     );
