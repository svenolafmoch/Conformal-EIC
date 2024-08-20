*
********************************************************************************
*
* hep-ph/0404111
* The Three-Loop Splitting Functions in QCD: The Singlet Case
* A. Vogt, S. Moch and J. Vermaseren
*
* This file contains the Mellin N results for singlet anomalous dimensions
* up to three loops as given in eqs.(3.5)--(3.13) and the results for the
* singlet splitting function in Bjorken x-space as given
* in eqs.(4.6)--(4.15). In Bjorken x-space all divergencies for x -> 1 are
* understood in the sense of +-distributions.
*
* The results for the quark-quark anomalous dimensions (splitting functions) as 
* given in eq.(2.4) are obtained by adding pure-singlet terms to the non-singlet 
* quantities g_ns^+ (P_ns^+). The reference for the latter quantities to three 
* loops is
* S. Moch, J. Vermaseren and A. Vogt,
* hep-ph/0403192 (to appear in Nucl. Phys. B) 
*
*
* Notation:
* S(R(m1,...,mk),N) denote harmonic sums
* H(R(m1,...,mk),x) denote harmonic polylogarithms
* z2 = Zeta(2), ... , z5=Zeta(5)
* den(x) = 1/x denotes denominators
*
********************************************************************************
*

L   gqg0 =
       + 4*theta( - 3 + N)*den(1 + N)*nf
       - 4*theta( - 3 + N)*den(2 + N)*nf
       - 2*theta( - 3 + N)*den(N)*nf
       - 2/3*delta( - 2 + N)*nf
      ;


L   ggg0 =
       - 11/3*theta( - 3 + N)*ca
       + 2/3*theta( - 3 + N)*nf
       - 4*theta( - 3 + N)*den( - 1 + N)*ca
       - 4*theta( - 3 + N)*den(1 + N)*ca
       + 4*theta( - 3 + N)*den(2 + N)*ca
       + 4*theta( - 3 + N)*den(N)*ca
       + 4*theta( - 3 + N)*S(R(1),N)*ca
       + 2/3*delta( - 2 + N)*nf
      ;

L   ggq0 =
       - 4*theta( - 3 + N)*den( - 1 + N)*cf
       - 2*theta( - 3 + N)*den(1 + N)*cf
       + 4*theta( - 3 + N)*den(N)*cf
       - 8/3*delta( - 2 + N)*cf
      ;


L   gqg1 =
       - 80/9*theta( - 3 + N)*den( - 1 + N)*ca*nf
       + 58*theta( - 3 + N)*den(1 + N)*cf*nf
       - 100*theta( - 3 + N)*den(1 + N)*ca*nf
       - 8*theta( - 3 + N)*den(1 + N)^2*cf*nf
       + 32*theta( - 3 + N)*den(1 + N)^2*ca*nf
       + 8*theta( - 3 + N)*den(1 + N)^3*cf*nf
       + 16*theta( - 3 + N)*den(1 + N)^3*ca*nf
       - 16*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*ca*nf
       + 16*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*ca*nf
       + 16*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*ca*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*nf
       - 40*theta( - 3 + N)*den(2 + N)*cf*nf
       + 872/9*theta( - 3 + N)*den(2 + N)*ca*nf
       + 16*theta( - 3 + N)*den(2 + N)^2*cf*nf
       + 176/3*theta( - 3 + N)*den(2 + N)^2*ca*nf
       - 16*theta( - 3 + N)*den(2 + N)^3*cf*nf
       + 16*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*nf
       + 16*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*ca*nf
       - 16*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*nf
       + 16*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*ca*nf
       - 16*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*nf
       + 16*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*ca*nf
       + 16*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf*nf
       - 28*theta( - 3 + N)*den(N)*cf*nf
       + 8*theta( - 3 + N)*den(N)*ca*nf
       + 6*theta( - 3 + N)*den(N)^2*cf*nf
       + 4*theta( - 3 + N)*den(N)^2*ca*nf
       - 4*theta( - 3 + N)*den(N)^3*cf*nf
       + 8*theta( - 3 + N)*den(N)^3*ca*nf
       + 8*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*nf
       + 8*theta( - 3 + N)*den(N)*S(R(-2),N)*ca*nf
       - 8*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*nf
       + 8*theta( - 3 + N)*den(N)*S(R(1,1),N)*ca*nf
       + 8*theta( - 3 + N)*den(N)*S(R(2),N)*cf*nf
       - 74/27*delta( - 2 + N)*cf*nf
       - 35/27*delta( - 2 + N)*ca*nf
      ;

L   ggg1 =
       + 2*theta( - 3 + N)*cf*nf
       + 8/3*theta( - 3 + N)*ca*nf
       - 32/3*theta( - 3 + N)*ca^2
       - 8/3*theta( - 3 + N)*den( - 1 + N)*cf*nf
       + 92/9*theta( - 3 + N)*den( - 1 + N)*ca*nf
       + 16*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*ca^2
       + 16*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*ca^2
       + 16*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*ca^2
       - 16*theta( - 3 + N)*den(1 + N)*cf*nf
       + 76/9*theta( - 3 + N)*den(1 + N)*ca*nf
       + 218/9*theta( - 3 + N)*den(1 + N)*ca^2
       - 20*theta( - 3 + N)*den(1 + N)^2*cf*nf
       - 8/3*theta( - 3 + N)*den(1 + N)^2*ca*nf
       + 44/3*theta( - 3 + N)*den(1 + N)^2*ca^2
       + 8*theta( - 3 + N)*den(1 + N)^3*cf*nf
       - 32*theta( - 3 + N)*den(1 + N)^3*ca^2
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*ca^2
       + 16*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*ca^2
       + 16*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*ca^2
       - 40/3*theta( - 3 + N)*den(2 + N)*cf*nf
       - 92/9*theta( - 3 + N)*den(2 + N)*ca*nf
       - 176/3*theta( - 3 + N)*den(2 + N)^2*ca^2
       + 16*theta( - 3 + N)*den(2 + N)^3*ca^2
       - 16*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*ca^2
       - 16*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*ca^2
       - 16*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*ca^2
       + 32*theta( - 3 + N)*den(N)*cf*nf
       - 76/9*theta( - 3 + N)*den(N)*ca*nf
       - 218/9*theta( - 3 + N)*den(N)*ca^2
       - 12*theta( - 3 + N)*den(N)^2*cf*nf
       - 8/3*theta( - 3 + N)*den(N)^2*ca*nf
       - 100/3*theta( - 3 + N)*den(N)^2*ca^2
       + 8*theta( - 3 + N)*den(N)^3*cf*nf
       - 16*theta( - 3 + N)*den(N)^3*ca^2
       - 16*theta( - 3 + N)*den(N)^2*S(R(1),N)*ca^2
       - 16*theta( - 3 + N)*den(N)*S(R(-2),N)*ca^2
       - 16*theta( - 3 + N)*den(N)*S(R(2),N)*ca^2
       + 8*theta( - 3 + N)*S(R(-3),N)*ca^2
       - 40/9*theta( - 3 + N)*S(R(1),N)*ca*nf
       + 268/9*theta( - 3 + N)*S(R(1),N)*ca^2
       - 16*theta( - 3 + N)*S(R(1,-2),N)*ca^2
       - 16*theta( - 3 + N)*S(R(1,2),N)*ca^2
       - 16*theta( - 3 + N)*S(R(2,1),N)*ca^2
       + 8*theta( - 3 + N)*S(R(3),N)*ca^2
       + 74/27*delta( - 2 + N)*cf*nf
       + 35/27*delta( - 2 + N)*ca*nf
      ;

L   ggq1 =
       - 4*theta( - 3 + N)*den( - 1 + N)*cf*ca
       + 80/9*theta( - 3 + N)*den( - 1 + N)*cf*nf
       + 16*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*cf*ca
       + 16*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*cf*ca
       + 88/3*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*ca
       - 16/3*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*nf
       - 24*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf^2
       - 16*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*ca
       + 16*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf^2
       + 16*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf*ca
       - 148/9*theta( - 3 + N)*den(1 + N)*cf*ca
       + 64/9*theta( - 3 + N)*den(1 + N)*cf*nf
       + 14*theta( - 3 + N)*den(1 + N)*cf^2
       - 20*theta( - 3 + N)*den(1 + N)^2*cf*ca
       + 14*theta( - 3 + N)*den(1 + N)^2*cf^2
       - 8*theta( - 3 + N)*den(1 + N)^3*cf*ca
       - 4*theta( - 3 + N)*den(1 + N)^3*cf^2
       + 8*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*ca
       + 8*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf*ca
       + 68/3*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*ca
       - 8/3*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*nf
       - 20*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^2
       - 8*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*ca
       + 8*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf^2
       + 8*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*ca
       - 176/9*theta( - 3 + N)*den(2 + N)*cf*ca
       - 32/3*theta( - 3 + N)*den(2 + N)^2*cf*ca
       - 76/9*theta( - 3 + N)*den(N)*cf*ca
       - 80/9*theta( - 3 + N)*den(N)*cf*nf
       + 10*theta( - 3 + N)*den(N)*cf^2
       - 48*theta( - 3 + N)*den(N)^2*cf*ca
       + 8*theta( - 3 + N)*den(N)^2*cf^2
       - 16*theta( - 3 + N)*den(N)^3*cf*ca
       + 8*theta( - 3 + N)*den(N)^3*cf^2
       - 16*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*ca
       - 16*theta( - 3 + N)*den(N)*S(R(-2),N)*cf*ca
       - 88/3*theta( - 3 + N)*den(N)*S(R(1),N)*cf*ca
       + 16/3*theta( - 3 + N)*den(N)*S(R(1),N)*cf*nf
       + 24*theta( - 3 + N)*den(N)*S(R(1),N)*cf^2
       + 16*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*ca
       - 16*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf^2
       - 16*theta( - 3 + N)*den(N)*S(R(2),N)*cf*ca
       - 376/27*delta( - 2 + N)*cf*ca
       + 104/27*delta( - 2 + N)*cf*nf
       + 112/27*delta( - 2 + N)*cf^2
      ;

L   gqqps1 =
       - 80/9*theta( - 3 + N)*den( - 1 + N)*cf*nf
       - 24*theta( - 3 + N)*den(1 + N)*cf*nf
       + 20*theta( - 3 + N)*den(1 + N)^2*cf*nf
       + 8*theta( - 3 + N)*den(1 + N)^3*cf*nf
       + 224/9*theta( - 3 + N)*den(2 + N)*cf*nf
       + 32/3*theta( - 3 + N)*den(2 + N)^2*cf*nf
       + 8*theta( - 3 + N)*den(N)*cf*nf
       + 4*theta( - 3 + N)*den(N)^2*cf*nf
       + 8*theta( - 3 + N)*den(N)^3*cf*nf
       - 40/27*delta( - 2 + N)*cf*nf
      ;

L   ggg2 =
       + 241/18*theta( - 3 + N)*cf*ca*nf
       - 11/9*theta( - 3 + N)*cf*nf^2
       - theta( - 3 + N)*cf^2*nf
       - 29/18*theta( - 3 + N)*ca*nf^2
       + 233/18*theta( - 3 + N)*ca^2*nf
       - 79/2*theta( - 3 + N)*ca^3
       - 64*theta( - 3 + N)*den( - 1 + N)*cf*ca*nf*z3
       + 30662/81*theta( - 3 + N)*den( - 1 + N)*cf*ca*nf
       + 1232/81*theta( - 3 + N)*den( - 1 + N)*cf*nf^2
       + 44/3*theta( - 3 + N)*den( - 1 + N)*cf^2*nf
       - 472/81*theta( - 3 + N)*den( - 1 + N)*ca*nf^2
       + 64*theta( - 3 + N)*den( - 1 + N)*ca^2*nf*z3
       - 19264/81*theta( - 3 + N)*den( - 1 + N)*ca^2*nf
       - 146182/81*theta( - 3 + N)*den( - 1 + N)*ca^3
       - 1376/27*theta( - 3 + N)*den( - 1 + N)^2*cf*ca*nf
       + 1136/27*theta( - 3 + N)*den( - 1 + N)^2*ca^2*nf
       + 6320/27*theta( - 3 + N)*den( - 1 + N)^2*ca^3
       - 128*theta( - 3 + N)*den( - 1 + N)^3*S(R(1,1), - 1 + N)*ca^3
       + 128*theta( - 3 + N)*den( - 1 + N)^3*S(R(2), - 1 + N)*ca^3
       + 224*theta( - 3 + N)*den( - 1 + N)^2*S(R(-3), - 1 + N)*ca^3
       - 64/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(-2), - 1 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(-2), - 1 + N)*ca^2*nf
       + 176/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(-2), - 1 + N)*ca^3
       - 128*theta( - 3 + N)*den( - 1 + N)^2*S(R(-2,1), - 1 + N)*ca^3
       + 32/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*cf*ca*nf
       - 176/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*ca^2*nf
       + 1072/9*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*ca^3
       - 64*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,-2), - 1 + N)*ca^3
       - 32/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,1), - 1 + N)*cf*ca*nf
       - 128*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,2), - 1 + N)*ca^3
       + 64/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(2), - 1 + N)*cf*ca*nf
       - 352/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(2), - 1 + N)*ca^3
       - 128*theta( - 3 + N)*den( - 1 + N)^2*S(R(2,1), - 1 + N)*ca^3
       + 160*theta( - 3 + N)*den( - 1 + N)^2*S(R(3), - 1 + N)*ca^3
       + 256*theta( - 3 + N)*den( - 1 + N)*S(R(-4), - 1 + N)*ca^3
       + 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*cf^2*nf
       + 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*ca^2*nf
       - 1408/3*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*ca^3
       - 256*theta( - 3 + N)*den( - 1 + N)*S(R(-3,1), - 1 + N)*ca^3
       + 704/9*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*cf*ca*nf
       - 304/3*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*ca^2*nf
       - 680/3*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*ca^3
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(-2,-2), - 1 + N)*ca^3
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1), - 1 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1), - 1 + N)*ca^2*nf
       + 176*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1), - 1 + N)*ca^3
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(-2,2), - 1 + N)*ca^3
       + 2084/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*ca*nf
       - 256/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*nf^2
       - 124/3*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf^2*nf
       + 104/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*ca*nf^2
       - 820/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*ca^2*nf
       + 1652/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*ca^3
       - 384*theta( - 3 + N)*den( - 1 + N)*S(R(1,-3), - 1 + N)*ca^3
       + 256/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*cf*ca*nf
       - 128/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*cf^2*nf
       - 32*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*ca^2*nf
       - 176*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*ca^3
       + 256*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2,1), - 1 + N)*ca^3
       - 16/9*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*ca*nf
       + 32/9*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*nf^2
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf^2*nf
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf*ca*nf
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf^2*nf
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf^2*nf
       - 128*theta( - 3 + N)*den( - 1 + N)*S(R(1,3), - 1 + N)*ca^3
       + 136/9*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf*ca*nf
       - 64/9*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf*nf^2
       - 80/3*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf^2*nf
       - 176/3*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*ca^2*nf
       + 1072/9*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*ca^3
       - 128*theta( - 3 + N)*den( - 1 + N)*S(R(2,-2), - 1 + N)*ca^3
       - 128*theta( - 3 + N)*den( - 1 + N)*S(R(2,2), - 1 + N)*ca^3
       + 128/3*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf*ca*nf
       - 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf^2*nf
       - 880/3*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*ca^3
       - 128*theta( - 3 + N)*den( - 1 + N)*S(R(3,1), - 1 + N)*ca^3
       + 128*theta( - 3 + N)*den( - 1 + N)*S(R(4), - 1 + N)*ca^3
       + 208*theta( - 3 + N)*den(1 + N)*cf*ca*nf*z3
       - 182/27*theta( - 3 + N)*den(1 + N)*cf*ca*nf
       - 64/9*theta( - 3 + N)*den(1 + N)*cf*nf^2
       - 144*theta( - 3 + N)*den(1 + N)*cf^2*nf*z3
       - 66*theta( - 3 + N)*den(1 + N)*cf^2*nf
       - 14/3*theta( - 3 + N)*den(1 + N)*ca*nf^2
       - 64*theta( - 3 + N)*den(1 + N)*ca^2*nf*z3
       - 7358/27*theta( - 3 + N)*den(1 + N)*ca^2*nf
       - 49678/27*theta( - 3 + N)*den(1 + N)*ca^3
       + 192*theta( - 3 + N)*den(1 + N)^2*cf*ca*nf*z3
       + 13804/27*theta( - 3 + N)*den(1 + N)^2*cf*ca*nf
       + 1088/27*theta( - 3 + N)*den(1 + N)^2*cf*nf^2
       - 96*theta( - 3 + N)*den(1 + N)^2*cf^2*nf*z3
       - 1082/3*theta( - 3 + N)*den(1 + N)^2*cf^2*nf
       + 86/27*theta( - 3 + N)*den(1 + N)^2*ca*nf^2
       - 96*theta( - 3 + N)*den(1 + N)^2*ca^2*nf*z3
       - 3958/27*theta( - 3 + N)*den(1 + N)^2*ca^2*nf
       + 9374/27*theta( - 3 + N)*den(1 + N)^2*ca^3
       - 4832/9*theta( - 3 + N)*den(1 + N)^3*cf*ca*nf
       - 232/9*theta( - 3 + N)*den(1 + N)^3*cf*nf^2
       + 396*theta( - 3 + N)*den(1 + N)^3*cf^2*nf
       + 16/9*theta( - 3 + N)*den(1 + N)^3*ca*nf^2
       + 2420/9*theta( - 3 + N)*den(1 + N)^3*ca^2*nf
       - 1052/3*theta( - 3 + N)*den(1 + N)^3*ca^3
       - 560/3*theta( - 3 + N)*den(1 + N)^4*cf*ca*nf
       - 16/3*theta( - 3 + N)*den(1 + N)^4*cf*nf^2
       - 88*theta( - 3 + N)*den(1 + N)^4*cf^2*nf
       - 152/3*theta( - 3 + N)*den(1 + N)^4*ca^2*nf
       + 464/3*theta( - 3 + N)*den(1 + N)^4*ca^3
       + 192*theta( - 3 + N)*den(1 + N)^5*cf*ca*nf
       + 32*theta( - 3 + N)*den(1 + N)^5*cf^2*nf
       - 16*theta( - 3 + N)*den(1 + N)^5*ca^2*nf
       - 448*theta( - 3 + N)*den(1 + N)^5*ca^3
       - 96*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf^2*nf
       + 24*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*ca^2*nf
       + 416*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*ca^3
       - 16*theta( - 3 + N)*den(1 + N)^3*S(R(-2),1 + N)*ca^2*nf
       + 64*theta( - 3 + N)*den(1 + N)^3*S(R(-2),1 + N)*ca^3
       - 152/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf*nf^2
       + 80*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf^2*nf
       + 212/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*ca^2*nf
       + 256/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*ca^3
       + 16*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf*ca*nf
       - 16*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf^2*nf
       - 128*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*ca^3
       - 64*theta( - 3 + N)*den(1 + N)^3*S(R(2),1 + N)*cf*ca*nf
       + 384*theta( - 3 + N)*den(1 + N)^3*S(R(2),1 + N)*ca^3
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*cf^2*nf
       - 48*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*ca^2*nf
       + 416*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*ca^3
       + 112*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf*ca*nf
       + 24*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*ca^2*nf
       - 160*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*ca^3
       - 128*theta( - 3 + N)*den(1 + N)^2*S(R(-2,1),1 + N)*ca^3
       - 596/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*ca*nf
       - 64/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*nf^2
       + 20*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf^2*nf
       - 32/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*ca*nf^2
       + 452/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*ca^2*nf
       - 128/3*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*ca^3
       - 256*theta( - 3 + N)*den(1 + N)^2*S(R(1,-2),1 + N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(1 + N)^2*S(R(1,-2),1 + N)*cf^2*nf
       + 96*theta( - 3 + N)*den(1 + N)^2*S(R(1,-2),1 + N)*ca^2*nf
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(1,-2),1 + N)*ca^3
       + 136/3*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf*ca*nf
       - 16/3*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf*nf^2
       - 64*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf^2*nf
       - 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf^2*nf
       + 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*cf^2*nf
       - 128*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*ca^3
       - 8/3*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf*nf^2
       + 128*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf^2*nf
       + 32/3*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*ca^2*nf
       - 464/3*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*ca^3
       - 128*theta( - 3 + N)*den(1 + N)^2*S(R(2,1),1 + N)*ca^3
       - 112*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf*ca*nf
       + 48*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf^2*nf
       + 24*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*ca^2*nf
       + 320*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*ca^3
       + 256*theta( - 3 + N)*den(1 + N)*S(R(-4),1 + N)*ca^3
       + 144*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf^2*nf
       - 72*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*ca^2*nf
       - 608*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*ca^3
       - 256*theta( - 3 + N)*den(1 + N)*S(R(-3,1),1 + N)*ca^3
       - 1520/3*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf*ca*nf
       + 1072/3*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf^2*nf
       + 1528/9*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*ca^2*nf
       - 3016/9*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*ca^3
       - 64*theta( - 3 + N)*den(1 + N)*S(R(-2,-2),1 + N)*ca^3
       - 64*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*ca^2*nf
       + 400*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*ca^3
       - 64*theta( - 3 + N)*den(1 + N)*S(R(-2,2),1 + N)*ca^3
       - 3428/9*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*ca*nf
       - 56/3*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*nf^2
       + 428/3*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^2*nf
       + 22/9*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*ca*nf^2
       + 1354/9*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*ca^2*nf
       + 124/9*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*ca^3
       - 384*theta( - 3 + N)*den(1 + N)*S(R(1,-3),1 + N)*ca^3
       - 288*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf^2*nf
       + 112*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*ca^2*nf
       + 48*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*ca^3
       + 256*theta( - 3 + N)*den(1 + N)*S(R(1,-2,1),1 + N)*ca^3
       + 12*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*ca*nf
       - 8/3*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*nf^2
       - 28/3*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf^2*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf^2*nf
       + 16*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf*ca*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf^2*nf
       - 128*theta( - 3 + N)*den(1 + N)*S(R(1,3),1 + N)*ca^3
       + 112/3*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*ca*nf
       + 16/3*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*nf^2
       + 40/3*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf^2*nf
       - 464/9*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*ca^2*nf
       + 200/9*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*ca^3
       - 128*theta( - 3 + N)*den(1 + N)*S(R(2,-2),1 + N)*ca^3
       - 128*theta( - 3 + N)*den(1 + N)*S(R(2,2),1 + N)*ca^3
       - 104*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf*ca*nf
       + 56*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf^2*nf
       + 100/3*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*ca^2*nf
       - 688/3*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*ca^3
       - 128*theta( - 3 + N)*den(1 + N)*S(R(3,1),1 + N)*ca^3
       + 128*theta( - 3 + N)*den(1 + N)*S(R(4),1 + N)*ca^3
       + 64*theta( - 3 + N)*den(2 + N)*cf*ca*nf*z3
       + 24526/81*theta( - 3 + N)*den(2 + N)*cf*ca*nf
       + 1360/81*theta( - 3 + N)*den(2 + N)*cf*nf^2
       - 232*theta( - 3 + N)*den(2 + N)*cf^2*nf
       + 472/81*theta( - 3 + N)*den(2 + N)*ca*nf^2
       - 64*theta( - 3 + N)*den(2 + N)*ca^2*nf*z3
       + 19264/81*theta( - 3 + N)*den(2 + N)*ca^2*nf
       + 146182/81*theta( - 3 + N)*den(2 + N)*ca^3
       + 112*theta( - 3 + N)*den(2 + N)^2*cf*ca*nf
       + 176/27*theta( - 3 + N)*den(2 + N)^2*cf*nf^2
       + 72*theta( - 3 + N)*den(2 + N)^2*cf^2*nf
       + 104/27*theta( - 3 + N)*den(2 + N)^2*ca*nf^2
       + 2612/9*theta( - 3 + N)*den(2 + N)^2*ca^2*nf
       + 10604/9*theta( - 3 + N)*den(2 + N)^2*ca^3
       - 1016/9*theta( - 3 + N)*den(2 + N)^3*cf*ca*nf
       - 64/9*theta( - 3 + N)*den(2 + N)^3*cf*nf^2
       - 224/3*theta( - 3 + N)*den(2 + N)^3*cf^2*nf
       - 496/9*theta( - 3 + N)*den(2 + N)^3*ca^2*nf
       + 12136/9*theta( - 3 + N)*den(2 + N)^3*ca^3
       + 32/3*theta( - 3 + N)*den(2 + N)^4*cf^2*nf
       - 704*theta( - 3 + N)*den(2 + N)^4*ca^3
       + 128*theta( - 3 + N)*den(2 + N)^5*ca^3
       - 256*theta( - 3 + N)*den(2 + N)^4*S(R(1),2 + N)*ca^3
       - 128*theta( - 3 + N)*den(2 + N)^3*S(R(-2),2 + N)*ca^3
       - 32*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*ca^2*nf
       + 1936/3*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*ca^3
       + 128*theta( - 3 + N)*den(2 + N)^3*S(R(1,1),2 + N)*ca^3
       - 192*theta( - 3 + N)*den(2 + N)^3*S(R(2),2 + N)*ca^3
       - 352*theta( - 3 + N)*den(2 + N)^2*S(R(-3),2 + N)*ca^3
       + 64*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*cf*ca*nf
       - 128/3*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*cf^2*nf
       - 64/3*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*ca^2*nf
       + 352*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*ca^3
       + 256*theta( - 3 + N)*den(2 + N)^2*S(R(-2,1),2 + N)*ca^3
       + 160/9*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*ca*nf
       + 32/9*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*nf^2
       + 64/3*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf^2*nf
       + 176/3*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*ca^2*nf
       - 1072/9*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*ca^3
       + 192*theta( - 3 + N)*den(2 + N)^2*S(R(1,-2),2 + N)*ca^3
       + 32/3*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf*ca*nf
       - 64/3*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf^2*nf
       + 128*theta( - 3 + N)*den(2 + N)^2*S(R(1,2),2 + N)*ca^3
       - 32/3*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*cf^2*nf
       + 352*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*ca^3
       + 128*theta( - 3 + N)*den(2 + N)^2*S(R(2,1),2 + N)*ca^3
       - 160*theta( - 3 + N)*den(2 + N)^2*S(R(3),2 + N)*ca^3
       - 256*theta( - 3 + N)*den(2 + N)*S(R(-4),2 + N)*ca^3
       - 32/3*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*cf*ca*nf
       - 64/3*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*cf^2*nf
       - 32/3*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*ca^2*nf
       + 1408/3*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*ca^3
       + 256*theta( - 3 + N)*den(2 + N)*S(R(-3,1),2 + N)*ca^3
       + 160/9*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*cf*ca*nf
       + 304/3*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*ca^2*nf
       + 680/3*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*ca^3
       + 64*theta( - 3 + N)*den(2 + N)*S(R(-2,-2),2 + N)*ca^3
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(-2,1),2 + N)*cf*ca*nf
       - 32/3*theta( - 3 + N)*den(2 + N)*S(R(-2,1),2 + N)*ca^2*nf
       - 176*theta( - 3 + N)*den(2 + N)*S(R(-2,1),2 + N)*ca^3
       + 64*theta( - 3 + N)*den(2 + N)*S(R(-2,2),2 + N)*ca^3
       + 1408/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*ca*nf
       - 176/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*nf^2
       - 72*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf^2*nf
       - 104/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*ca*nf^2
       + 820/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*ca^2*nf
       - 1652/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*ca^3
       + 384*theta( - 3 + N)*den(2 + N)*S(R(1,-3),2 + N)*ca^3
       - 256/3*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*cf*ca*nf
       + 128/3*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*cf^2*nf
       + 32*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*ca^2*nf
       + 176*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*ca^3
       - 256*theta( - 3 + N)*den(2 + N)*S(R(1,-2,1),2 + N)*ca^3
       + 16/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*ca*nf
       - 32/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*nf^2
       - 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf^2*nf
       - 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf^2*nf
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,2),2 + N)*cf*ca*nf
       - 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,2),2 + N)*cf^2*nf
       + 128*theta( - 3 + N)*den(2 + N)*S(R(1,3),2 + N)*ca^3
       + 296/9*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf*ca*nf
       + 64/9*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf*nf^2
       + 224/3*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf^2*nf
       + 176/3*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*ca^2*nf
       - 1072/9*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*ca^3
       + 128*theta( - 3 + N)*den(2 + N)*S(R(2,-2),2 + N)*ca^3
       + 128*theta( - 3 + N)*den(2 + N)*S(R(2,2),2 + N)*ca^3
       - 128/3*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*cf^2*nf
       + 880/3*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*ca^3
       + 128*theta( - 3 + N)*den(2 + N)*S(R(3,1),2 + N)*ca^3
       - 128*theta( - 3 + N)*den(2 + N)*S(R(4),2 + N)*ca^3
       - 208*theta( - 3 + N)*den(N)*cf*ca*nf*z3
       - 18214/27*theta( - 3 + N)*den(N)*cf*ca*nf
       - 224/9*theta( - 3 + N)*den(N)*cf*nf^2
       + 144*theta( - 3 + N)*den(N)*cf^2*nf*z3
       + 850/3*theta( - 3 + N)*den(N)*cf^2*nf
       + 14/3*theta( - 3 + N)*den(N)*ca*nf^2
       + 64*theta( - 3 + N)*den(N)*ca^2*nf*z3
       + 7358/27*theta( - 3 + N)*den(N)*ca^2*nf
       + 49678/27*theta( - 3 + N)*den(N)*ca^3
       + 192*theta( - 3 + N)*den(N)^2*cf*ca*nf*z3
       + 3220/27*theta( - 3 + N)*den(N)^2*cf*ca*nf
       - 352/27*theta( - 3 + N)*den(N)^2*cf*nf^2
       - 96*theta( - 3 + N)*den(N)^2*cf^2*nf*z3
       - 254*theta( - 3 + N)*den(N)^2*cf^2*nf
       + 152/27*theta( - 3 + N)*den(N)^2*ca*nf^2
       - 96*theta( - 3 + N)*den(N)^2*ca^2*nf*z3
       + 7400/27*theta( - 3 + N)*den(N)^2*ca^2*nf
       + 30866/27*theta( - 3 + N)*den(N)^2*ca^3
       - 4028/9*theta( - 3 + N)*den(N)^3*cf*ca*nf
       - 184/9*theta( - 3 + N)*den(N)^3*cf*nf^2
       + 328/3*theta( - 3 + N)*den(N)^3*cf^2*nf
       + 16/9*theta( - 3 + N)*den(N)^3*ca*nf^2
       + 176/3*theta( - 3 + N)*den(N)^3*ca^2*nf
       + 9476/9*theta( - 3 + N)*den(N)^3*ca^3
       + 232/3*theta( - 3 + N)*den(N)^4*cf*ca*nf
       - 16/3*theta( - 3 + N)*den(N)^4*cf*nf^2
       - 24*theta( - 3 + N)*den(N)^4*cf^2*nf
       + 32*theta( - 3 + N)*den(N)^4*ca^2*nf
       - 368*theta( - 3 + N)*den(N)^4*ca^3
       - 128*theta( - 3 + N)*den(N)^5*cf*ca*nf
       + 32*theta( - 3 + N)*den(N)^5*cf^2*nf
       + 128*theta( - 3 + N)*den(N)^5*ca^3
       - 96*theta( - 3 + N)*den(N)^4*S(R(1),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)^4*S(R(1),N)*cf^2*nf
       + 24*theta( - 3 + N)*den(N)^4*S(R(1),N)*ca^2*nf
       + 160*theta( - 3 + N)*den(N)^4*S(R(1),N)*ca^3
       + 128*theta( - 3 + N)*den(N)^3*S(R(-2),N)*cf*ca*nf
       - 128*theta( - 3 + N)*den(N)^3*S(R(-2),N)*cf^2*nf
       - 80*theta( - 3 + N)*den(N)^3*S(R(-2),N)*ca^2*nf
       + 192*theta( - 3 + N)*den(N)^3*S(R(-2),N)*ca^3
       + 64/3*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf*nf^2
       + 64*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf^2*nf
       + 4*theta( - 3 + N)*den(N)^3*S(R(1),N)*ca^2*nf
       + 464*theta( - 3 + N)*den(N)^3*S(R(1),N)*ca^3
       + 16*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf*ca*nf
       - 16*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf^2*nf
       + 128*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*ca^3
       - 64*theta( - 3 + N)*den(N)^3*S(R(2),N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(N)^3*S(R(2),N)*ca^3
       + 64*theta( - 3 + N)*den(N)^2*S(R(-3),N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(N)^2*S(R(-3),N)*cf^2*nf
       - 48*theta( - 3 + N)*den(N)^2*S(R(-3),N)*ca^2*nf
       - 160*theta( - 3 + N)*den(N)^2*S(R(-3),N)*ca^3
       - 208*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf^2*nf
       + 136*theta( - 3 + N)*den(N)^2*S(R(-2),N)*ca^2*nf
       + 272*theta( - 3 + N)*den(N)^2*S(R(-2),N)*ca^3
       + 256*theta( - 3 + N)*den(N)^2*S(R(-2,1),N)*ca^3
       - 2696/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*ca*nf
       - 16/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*nf^2
       + 64*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf^2*nf
       - 32/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*ca*nf^2
       + 460/3*theta( - 3 + N)*den(N)^2*S(R(1),N)*ca^2*nf
       - 784/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*ca^3
       - 256*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*cf^2*nf
       + 96*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*ca^2*nf
       + 320*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*ca^3
       + 88/3*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf*ca*nf
       - 16/3*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf*nf^2
       - 48*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf^2*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf^2*nf
       + 32*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*cf^2*nf
       + 128*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*ca^3
       + 16/3*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf*nf^2
       + 80*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf^2*nf
       + 32/3*theta( - 3 + N)*den(N)^2*S(R(2),N)*ca^2*nf
       + 688/3*theta( - 3 + N)*den(N)^2*S(R(2),N)*ca^3
       + 128*theta( - 3 + N)*den(N)^2*S(R(2,1),N)*ca^3
       - 112*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf*ca*nf
       + 48*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf^2*nf
       + 24*theta( - 3 + N)*den(N)^2*S(R(3),N)*ca^2*nf
       - 256*theta( - 3 + N)*den(N)*S(R(-4),N)*ca^3
       - 144*theta( - 3 + N)*den(N)*S(R(-3),N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(N)*S(R(-3),N)*cf^2*nf
       + 72*theta( - 3 + N)*den(N)*S(R(-3),N)*ca^2*nf
       + 608*theta( - 3 + N)*den(N)*S(R(-3),N)*ca^3
       + 256*theta( - 3 + N)*den(N)*S(R(-3,1),N)*ca^3
       + 1232/3*theta( - 3 + N)*den(N)*S(R(-2),N)*cf*ca*nf
       - 1072/3*theta( - 3 + N)*den(N)*S(R(-2),N)*cf^2*nf
       - 1528/9*theta( - 3 + N)*den(N)*S(R(-2),N)*ca^2*nf
       + 3016/9*theta( - 3 + N)*den(N)*S(R(-2),N)*ca^3
       + 64*theta( - 3 + N)*den(N)*S(R(-2,-2),N)*ca^3
       + 64*theta( - 3 + N)*den(N)*S(R(-2,1),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)*S(R(-2,1),N)*ca^2*nf
       - 400*theta( - 3 + N)*den(N)*S(R(-2,1),N)*ca^3
       + 64*theta( - 3 + N)*den(N)*S(R(-2,2),N)*ca^3
       + 2264/9*theta( - 3 + N)*den(N)*S(R(1),N)*cf*ca*nf
       + 104/3*theta( - 3 + N)*den(N)*S(R(1),N)*cf*nf^2
       - 88/3*theta( - 3 + N)*den(N)*S(R(1),N)*cf^2*nf
       - 22/9*theta( - 3 + N)*den(N)*S(R(1),N)*ca*nf^2
       - 1354/9*theta( - 3 + N)*den(N)*S(R(1),N)*ca^2*nf
       - 124/9*theta( - 3 + N)*den(N)*S(R(1),N)*ca^3
       + 384*theta( - 3 + N)*den(N)*S(R(1,-3),N)*ca^3
       + 288*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf*ca*nf
       - 128*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf^2*nf
       - 112*theta( - 3 + N)*den(N)*S(R(1,-2),N)*ca^2*nf
       - 48*theta( - 3 + N)*den(N)*S(R(1,-2),N)*ca^3
       - 256*theta( - 3 + N)*den(N)*S(R(1,-2,1),N)*ca^3
       - 12*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*ca*nf
       + 8/3*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*nf^2
       + 28/3*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf^2*nf
       + 16*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf*ca*nf
       - 16*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf^2*nf
       - 16*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf^2*nf
       + 128*theta( - 3 + N)*den(N)*S(R(1,3),N)*ca^3
       - 256/3*theta( - 3 + N)*den(N)*S(R(2),N)*cf*ca*nf
       - 16/3*theta( - 3 + N)*den(N)*S(R(2),N)*cf*nf^2
       - 184/3*theta( - 3 + N)*den(N)*S(R(2),N)*cf^2*nf
       + 464/9*theta( - 3 + N)*den(N)*S(R(2),N)*ca^2*nf
       - 200/9*theta( - 3 + N)*den(N)*S(R(2),N)*ca^3
       + 128*theta( - 3 + N)*den(N)*S(R(2,-2),N)*ca^3
       + 128*theta( - 3 + N)*den(N)*S(R(2,2),N)*ca^3
       + 104*theta( - 3 + N)*den(N)*S(R(3),N)*cf*ca*nf
       - 56*theta( - 3 + N)*den(N)*S(R(3),N)*cf^2*nf
       - 100/3*theta( - 3 + N)*den(N)*S(R(3),N)*ca^2*nf
       + 688/3*theta( - 3 + N)*den(N)*S(R(3),N)*ca^3
       + 128*theta( - 3 + N)*den(N)*S(R(3,1),N)*ca^3
       - 128*theta( - 3 + N)*den(N)*S(R(4),N)*ca^3
       + 64*theta( - 3 + N)*S(R(-5),N)*ca^3
       - 128*theta( - 3 + N)*S(R(-4,1),N)*ca^3
       - 160/9*theta( - 3 + N)*S(R(-3),N)*ca^2*nf
       + 1072/9*theta( - 3 + N)*S(R(-3),N)*ca^3
       - 64*theta( - 3 + N)*S(R(-3,-2),N)*ca^3
       - 32*theta( - 3 + N)*S(R(-3,2),N)*ca^3
       - 64*theta( - 3 + N)*S(R(-2,-3),N)*ca^3
       + 32/3*theta( - 3 + N)*S(R(-2,-2),N)*ca^2*nf
       - 176/3*theta( - 3 + N)*S(R(-2,-2),N)*ca^3
       + 64*theta( - 3 + N)*S(R(-2,-2,1),N)*ca^3
       + 64*theta( - 3 + N)*S(R(-2,1,-2),N)*ca^3
       + 32*theta( - 3 + N)*S(R(1),N)*cf*ca*nf*z3
       - 110/3*theta( - 3 + N)*S(R(1),N)*cf*ca*nf
       - 16/27*theta( - 3 + N)*S(R(1),N)*ca*nf^2
       - 32*theta( - 3 + N)*S(R(1),N)*ca^2*nf*z3
       - 836/27*theta( - 3 + N)*S(R(1),N)*ca^2*nf
       + 490/3*theta( - 3 + N)*S(R(1),N)*ca^3
       - 256*theta( - 3 + N)*S(R(1,-4),N)*ca^3
       + 256*theta( - 3 + N)*S(R(1,-3,1),N)*ca^3
       + 320/9*theta( - 3 + N)*S(R(1,-2),N)*ca^2*nf
       - 2144/9*theta( - 3 + N)*S(R(1,-2),N)*ca^3
       + 64*theta( - 3 + N)*S(R(1,-2,-2),N)*ca^3
       + 64*theta( - 3 + N)*S(R(1,-2,2),N)*ca^3
       + 384*theta( - 3 + N)*S(R(1,1,-3),N)*ca^3
       - 256*theta( - 3 + N)*S(R(1,1,-2,1),N)*ca^3
       + 128*theta( - 3 + N)*S(R(1,1,3),N)*ca^3
       + 320/9*theta( - 3 + N)*S(R(1,2),N)*ca^2*nf
       - 2144/9*theta( - 3 + N)*S(R(1,2),N)*ca^3
       + 128*theta( - 3 + N)*S(R(1,2,-2),N)*ca^3
       + 128*theta( - 3 + N)*S(R(1,2,2),N)*ca^3
       - 16/3*theta( - 3 + N)*S(R(1,3),N)*ca^2*nf
       + 88/3*theta( - 3 + N)*S(R(1,3),N)*ca^3
       + 128*theta( - 3 + N)*S(R(1,3,1),N)*ca^3
       - 128*theta( - 3 + N)*S(R(1,4),N)*ca^3
       + 8/3*theta( - 3 + N)*S(R(2),N)*ca^2*nf
       - 8/3*theta( - 3 + N)*S(R(2),N)*ca^3
       - 288*theta( - 3 + N)*S(R(2,-3),N)*ca^3
       + 192*theta( - 3 + N)*S(R(2,-2,1),N)*ca^3
       + 320/9*theta( - 3 + N)*S(R(2,1),N)*ca^2*nf
       - 2144/9*theta( - 3 + N)*S(R(2,1),N)*ca^3
       + 128*theta( - 3 + N)*S(R(2,1,-2),N)*ca^3
       + 128*theta( - 3 + N)*S(R(2,1,2),N)*ca^3
       + 128*theta( - 3 + N)*S(R(2,2,1),N)*ca^3
       - 160*theta( - 3 + N)*S(R(2,3),N)*ca^3
       - 160/9*theta( - 3 + N)*S(R(3),N)*ca^2*nf
       + 1072/9*theta( - 3 + N)*S(R(3),N)*ca^3
       - 64*theta( - 3 + N)*S(R(3,-2),N)*ca^3
       - 16/3*theta( - 3 + N)*S(R(3,1),N)*ca^2*nf
       + 88/3*theta( - 3 + N)*S(R(3,1),N)*ca^3
       + 128*theta( - 3 + N)*S(R(3,1,1),N)*ca^3
       - 160*theta( - 3 + N)*S(R(3,2),N)*ca^3
       - 128*theta( - 3 + N)*S(R(4,1),N)*ca^3
       + 64*theta( - 3 + N)*S(R(5),N)*ca^3
       + 104/3*delta( - 2 + N)*cf*ca*nf*z3
       - 139/9*delta( - 2 + N)*cf*ca*nf
       + 173/243*delta( - 2 + N)*cf*nf^2
       - 32/3*delta( - 2 + N)*cf^2*nf*z3
       + 2155/243*delta( - 2 + N)*cf^2*nf
       - 1058/243*delta( - 2 + N)*ca*nf^2
       - 24*delta( - 2 + N)*ca^2*nf*z3
       + 3589/162*delta( - 2 + N)*ca^2*nf
      ;

L   gqg2 =
       + 64*theta( - 3 + N)*den( - 1 + N)*cf*ca*nf*z3
       - 220/3*theta( - 3 + N)*den( - 1 + N)*cf*ca*nf
       - 1232/81*theta( - 3 + N)*den( - 1 + N)*cf*nf^2
       + 424/81*theta( - 3 + N)*den( - 1 + N)*ca*nf^2
       - 64*theta( - 3 + N)*den( - 1 + N)*ca^2*nf*z3
       + 9404/27*theta( - 3 + N)*den( - 1 + N)*ca^2*nf
       - 896/27*theta( - 3 + N)*den( - 1 + N)^2*ca^2*nf
       + 320/9*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*ca^2*nf
       + 32/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,1), - 1 + N)*ca^2*nf
       + 128/3*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*ca^2*nf
       + 832/9*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*ca^2*nf
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1), - 1 + N)*ca^2*nf
       + 1484/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*ca*nf
       - 644/9*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*ca^2*nf
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*ca^2*nf
       - 368/9*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*ca*nf
       + 160/9*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*ca^2*nf
       + 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf*ca*nf
       - 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*ca^2*nf
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*ca^2*nf
       + 368/9*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf*ca*nf
       + 24*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*ca^2*nf
       - 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf*ca*nf
       + 32*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*ca^2*nf
       - 96*theta( - 3 + N)*den(1 + N)*cf*ca*nf*z3
       - 2303/2*theta( - 3 + N)*den(1 + N)*cf*ca*nf
       + 89/27*theta( - 3 + N)*den(1 + N)*cf*nf^2
       + 288*theta( - 3 + N)*den(1 + N)*cf^2*nf*z3
       + 641/2*theta( - 3 + N)*den(1 + N)*cf^2*nf
       + 1348/27*theta( - 3 + N)*den(1 + N)*ca*nf^2
       - 192*theta( - 3 + N)*den(1 + N)*ca^2*nf*z3
       + 102608/27*theta( - 3 + N)*den(1 + N)*ca^2*nf
       - 480*theta( - 3 + N)*den(1 + N)^2*cf*ca*nf*z3
       + 24569/27*theta( - 3 + N)*den(1 + N)^2*cf*ca*nf
       - 5270/27*theta( - 3 + N)*den(1 + N)^2*cf*nf^2
       + 96*theta( - 3 + N)*den(1 + N)^2*cf^2*nf*z3
       + 57*theta( - 3 + N)*den(1 + N)^2*cf^2*nf
       + 560/27*theta( - 3 + N)*den(1 + N)^2*ca*nf^2
       + 384*theta( - 3 + N)*den(1 + N)^2*ca^2*nf*z3
       - 30962/27*theta( - 3 + N)*den(1 + N)^2*ca^2*nf
       + 3190/9*theta( - 3 + N)*den(1 + N)^3*cf*ca*nf
       + 884/9*theta( - 3 + N)*den(1 + N)^3*cf*nf^2
       - 126*theta( - 3 + N)*den(1 + N)^3*cf^2*nf
       - 104/3*theta( - 3 + N)*den(1 + N)^3*ca*nf^2
       - 200*theta( - 3 + N)*den(1 + N)^3*ca^2*nf
       + 344/3*theta( - 3 + N)*den(1 + N)^4*cf*ca*nf
       + 160/3*theta( - 3 + N)*den(1 + N)^4*cf*nf^2
       + 64*theta( - 3 + N)*den(1 + N)^4*cf^2*nf
       + 64/3*theta( - 3 + N)*den(1 + N)^4*ca*nf^2
       + 56/3*theta( - 3 + N)*den(1 + N)^4*ca^2*nf
       + 64*theta( - 3 + N)*den(1 + N)^5*cf*ca*nf
       - 64*theta( - 3 + N)*den(1 + N)^5*cf*nf^2
       + 240*theta( - 3 + N)*den(1 + N)^5*ca^2*nf
       - 128*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf*ca*nf
       - 112*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf^2*nf
       - 72*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*ca^2*nf
       - 32*theta( - 3 + N)*den(1 + N)^3*S(R(-2),1 + N)*cf*ca*nf
       - 16*theta( - 3 + N)*den(1 + N)^3*S(R(-2),1 + N)*ca^2*nf
       - 296*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf*ca*nf
       + 160*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf^2*nf
       - 108*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*ca^2*nf
       + 96*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf*ca*nf
       + 160*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf^2*nf
       - 128*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(1 + N)^3*S(R(2),1 + N)*cf*ca*nf
       - 128*theta( - 3 + N)*den(1 + N)^3*S(R(2),1 + N)*cf^2*nf
       - 96*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*cf^2*nf
       - 208*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*ca^2*nf
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf^2*nf
       - 216*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*ca^2*nf
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(-2,1),1 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(-2,1),1 + N)*ca^2*nf
       + 1100/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*ca*nf
       + 112/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*nf^2
       - 416*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf^2*nf
       + 8/3*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*ca*nf^2
       + 1256/3*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*ca^2*nf
       + 128*theta( - 3 + N)*den(1 + N)^2*S(R(1,-2),1 + N)*cf*ca*nf
       - 160*theta( - 3 + N)*den(1 + N)^2*S(R(1,-2),1 + N)*ca^2*nf
       + 496/3*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf*ca*nf
       - 16/3*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf*nf^2
       - 104*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf^2*nf
       - 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*ca^2*nf
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf*ca*nf
       - 80*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf^2*nf
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*ca^2*nf
       + 192*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*cf^2*nf
       - 128*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*ca^2*nf
       - 72*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf*ca*nf
       + 80*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf^2*nf
       - 96*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*ca^2*nf
       + 128*theta( - 3 + N)*den(1 + N)^2*S(R(2,1),1 + N)*cf^2*nf
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf*ca*nf
       - 96*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf^2*nf
       - 200*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*ca^2*nf
       - 176*theta( - 3 + N)*den(1 + N)*S(R(-4),1 + N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(1 + N)*S(R(-4),1 + N)*cf^2*nf
       - 48*theta( - 3 + N)*den(1 + N)*S(R(-4),1 + N)*ca^2*nf
       - 64*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf*ca*nf
       + 288*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf^2*nf
       - 32/3*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*ca*nf^2
       + 1016/3*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*ca^2*nf
       + 192*theta( - 3 + N)*den(1 + N)*S(R(-3,1),1 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(-3,1),1 + N)*ca^2*nf
       + 152*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf*ca*nf
       - 368*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf^2*nf
       + 112/9*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*ca*nf^2
       + 5864/9*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*ca^2*nf
       - 96*theta( - 3 + N)*den(1 + N)*S(R(-2,-2),1 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(-2,-2),1 + N)*cf^2*nf
       + 96*theta( - 3 + N)*den(1 + N)*S(R(-2,-2),1 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*cf*ca*nf
       - 128*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*ca^2*nf
       - 64*theta( - 3 + N)*den(1 + N)*S(R(-2,1,1),1 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(-2,1,1),1 + N)*ca^2*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(-2,2),1 + N)*cf*ca*nf
       + 96*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*ca*nf*z3
       - 12892/27*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*ca*nf
       + 148/27*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*nf^2
       + 282*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^2*nf
       - 112/27*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*ca*nf^2
       - 96*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*ca^2*nf*z3
       + 7750/27*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*ca^2*nf
       + 192*theta( - 3 + N)*den(1 + N)*S(R(1,-3),1 + N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(1 + N)*S(R(1,-3),1 + N)*cf^2*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,-3),1 + N)*ca^2*nf
       + 224*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf*ca*nf
       - 256*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf^2*nf
       + 80*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(1 + N)*S(R(1,-2,1),1 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(1 + N)*S(R(1,-2,1),1 + N)*ca^2*nf
       - 4568/9*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*ca*nf
       - 112/9*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*nf^2
       + 296*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf^2*nf
       + 112/9*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*ca*nf^2
       + 1904/9*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(1 + N)*S(R(1,1,-2),1 + N)*cf*ca*nf
       + 192*theta( - 3 + N)*den(1 + N)*S(R(1,1,-2),1 + N)*ca^2*nf
       - 424/3*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf*ca*nf
       + 16/3*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf*nf^2
       + 80*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf^2*nf
       - 16/3*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*ca*nf^2
       + 184/3*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*ca^2*nf
       - 128*theta( - 3 + N)*den(1 + N)*S(R(1,1,1,1),1 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,1,1,1),1 + N)*cf^2*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,1,1,1),1 + N)*ca^2*nf
       - 96*theta( - 3 + N)*den(1 + N)*S(R(1,1,2),1 + N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(1 + N)*S(R(1,1,2),1 + N)*cf^2*nf
       + 128*theta( - 3 + N)*den(1 + N)*S(R(1,1,2),1 + N)*ca^2*nf
       - 88*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf*ca*nf
       - 48*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf^2*nf
       + 16/3*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*ca*nf^2
       + 392/3*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*ca^2*nf
       + 96*theta( - 3 + N)*den(1 + N)*S(R(1,2,1),1 + N)*cf*ca*nf
       - 160*theta( - 3 + N)*den(1 + N)*S(R(1,2,1),1 + N)*cf^2*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,2,1),1 + N)*ca^2*nf
       + 32*theta( - 3 + N)*den(1 + N)*S(R(1,3),1 + N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(1 + N)*S(R(1,3),1 + N)*cf^2*nf
       - 32*theta( - 3 + N)*den(1 + N)*S(R(1,3),1 + N)*ca^2*nf
       + 1888/9*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*ca*nf
       + 128/9*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*nf^2
       - 352*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf^2*nf
       + 688/3*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*ca^2*nf
       + 96*theta( - 3 + N)*den(1 + N)*S(R(2,-2),1 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(2,-2),1 + N)*cf^2*nf
       - 32*theta( - 3 + N)*den(1 + N)*S(R(2,-2),1 + N)*ca^2*nf
       + 152*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*cf*ca*nf
       - 208*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*cf^2*nf
       - 16/3*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*ca*nf^2
       + 184/3*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*ca^2*nf
       + 192*theta( - 3 + N)*den(1 + N)*S(R(2,1,1),1 + N)*cf*ca*nf
       - 192*theta( - 3 + N)*den(1 + N)*S(R(2,1,1),1 + N)*cf^2*nf
       + 192*theta( - 3 + N)*den(1 + N)*S(R(2,2),1 + N)*cf^2*nf
       - 64*theta( - 3 + N)*den(1 + N)*S(R(2,2),1 + N)*ca^2*nf
       + 488/3*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf*ca*nf
       - 32/3*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf*nf^2
       + 48*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf^2*nf
       - 16/3*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*ca*nf^2
       + 148/3*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(1 + N)*S(R(3,1),1 + N)*cf*ca*nf
       + 224*theta( - 3 + N)*den(1 + N)*S(R(3,1),1 + N)*cf^2*nf
       + 96*theta( - 3 + N)*den(1 + N)*S(R(3,1),1 + N)*ca^2*nf
       + 48*theta( - 3 + N)*den(1 + N)*S(R(4),1 + N)*cf*ca*nf
       - 96*theta( - 3 + N)*den(1 + N)*S(R(4),1 + N)*cf^2*nf
       - 80*theta( - 3 + N)*den(1 + N)*S(R(4),1 + N)*ca^2*nf
       - 496*theta( - 3 + N)*den(2 + N)*cf*ca*nf*z3
       + 3997/3*theta( - 3 + N)*den(2 + N)*cf*ca*nf
       - 18424/81*theta( - 3 + N)*den(2 + N)*cf*nf^2
       - 81*theta( - 3 + N)*den(2 + N)*cf^2*nf
       - 2800/81*theta( - 3 + N)*den(2 + N)*ca*nf^2
       + 496*theta( - 3 + N)*den(2 + N)*ca^2*nf*z3
       - 103040/27*theta( - 3 + N)*den(2 + N)*ca^2*nf
       - 96*theta( - 3 + N)*den(2 + N)^2*cf*ca*nf*z3
       + 11332/27*theta( - 3 + N)*den(2 + N)^2*cf*ca*nf
       - 560/9*theta( - 3 + N)*den(2 + N)^2*cf*nf^2
       + 348*theta( - 3 + N)*den(2 + N)^2*cf^2*nf
       - 1832/27*theta( - 3 + N)*den(2 + N)^2*ca*nf^2
       + 96*theta( - 3 + N)*den(2 + N)^2*ca^2*nf*z3
       - 3040*theta( - 3 + N)*den(2 + N)^2*ca^2*nf
       + 2204/9*theta( - 3 + N)*den(2 + N)^3*cf*ca*nf
       + 352/9*theta( - 3 + N)*den(2 + N)^3*cf*nf^2
       - 352*theta( - 3 + N)*den(2 + N)^3*cf^2*nf
       - 64/9*theta( - 3 + N)*den(2 + N)^3*ca*nf^2
       - 10244/9*theta( - 3 + N)*den(2 + N)^3*ca^2*nf
       + 1264/3*theta( - 3 + N)*den(2 + N)^4*cf*ca*nf
       + 32*theta( - 3 + N)*den(2 + N)^4*cf*nf^2
       + 240*theta( - 3 + N)*den(2 + N)^4*cf^2*nf
       + 96*theta( - 3 + N)*den(2 + N)^4*ca^2*nf
       - 128*theta( - 3 + N)*den(2 + N)^5*cf^2*nf
       - 64*theta( - 3 + N)*den(2 + N)^4*S(R(1),2 + N)*cf*ca*nf
       + 224*theta( - 3 + N)*den(2 + N)^4*S(R(1),2 + N)*cf^2*nf
       + 96*theta( - 3 + N)*den(2 + N)^4*S(R(1),2 + N)*ca^2*nf
       + 128*theta( - 3 + N)*den(2 + N)^3*S(R(-2),2 + N)*cf^2*nf
       - 1384/3*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*cf*ca*nf
       - 176*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*cf^2*nf
       - 248/3*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*ca^2*nf
       + 64*theta( - 3 + N)*den(2 + N)^3*S(R(1,1),2 + N)*cf*ca*nf
       - 192*theta( - 3 + N)*den(2 + N)^3*S(R(1,1),2 + N)*cf^2*nf
       + 192*theta( - 3 + N)*den(2 + N)^3*S(R(2),2 + N)*cf^2*nf
       + 192*theta( - 3 + N)*den(2 + N)^2*S(R(-3),2 + N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(2 + N)^2*S(R(-3),2 + N)*cf^2*nf
       + 32*theta( - 3 + N)*den(2 + N)^2*S(R(-3),2 + N)*ca^2*nf
       + 80*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*cf*ca*nf
       - 128*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*cf^2*nf
       - 400*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(2 + N)^2*S(R(-2,1),2 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(2 + N)^2*S(R(-2,1),2 + N)*ca^2*nf
       - 392*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*ca*nf
       - 112/9*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*nf^2
       + 244*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf^2*nf
       + 620/9*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(2 + N)^2*S(R(1,-2),2 + N)*cf*ca*nf
       + 192*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf*ca*nf
       + 16/3*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf*nf^2
       + 96*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf^2*nf
       - 848/3*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*ca^2*nf
       - 64*theta( - 3 + N)*den(2 + N)^2*S(R(1,1,1),2 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(2 + N)^2*S(R(1,1,1),2 + N)*cf^2*nf
       - 96*theta( - 3 + N)*den(2 + N)^2*S(R(1,2),2 + N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(2 + N)^2*S(R(1,2),2 + N)*cf^2*nf
       - 392*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*cf*ca*nf
       - 80*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*cf^2*nf
       + 40/3*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*ca^2*nf
       + 32*theta( - 3 + N)*den(2 + N)^2*S(R(2,1),2 + N)*cf*ca*nf
       - 160*theta( - 3 + N)*den(2 + N)^2*S(R(2,1),2 + N)*cf^2*nf
       + 96*theta( - 3 + N)*den(2 + N)^2*S(R(3),2 + N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(2 + N)^2*S(R(3),2 + N)*cf^2*nf
       - 64*theta( - 3 + N)*den(2 + N)^2*S(R(3),2 + N)*ca^2*nf
       + 176*theta( - 3 + N)*den(2 + N)*S(R(-4),2 + N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(2 + N)*S(R(-4),2 + N)*cf^2*nf
       + 48*theta( - 3 + N)*den(2 + N)*S(R(-4),2 + N)*ca^2*nf
       + 48*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*cf*ca*nf
       - 160*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*cf^2*nf
       + 32/3*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*ca*nf^2
       - 976/3*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(2 + N)*S(R(-3,1),2 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(-3,1),2 + N)*ca^2*nf
       + 144*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*cf*ca*nf
       + 80*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*cf^2*nf
       - 112/9*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*ca*nf^2
       - 2216/3*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*ca^2*nf
       + 96*theta( - 3 + N)*den(2 + N)*S(R(-2,-2),2 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(-2,-2),2 + N)*cf^2*nf
       - 96*theta( - 3 + N)*den(2 + N)*S(R(-2,-2),2 + N)*ca^2*nf
       + 80*theta( - 3 + N)*den(2 + N)*S(R(-2,1),2 + N)*cf*ca*nf
       + 256/3*theta( - 3 + N)*den(2 + N)*S(R(-2,1),2 + N)*ca^2*nf
       + 64*theta( - 3 + N)*den(2 + N)*S(R(-2,1,1),2 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(-2,1,1),2 + N)*ca^2*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(-2,2),2 + N)*cf*ca*nf
       - 96*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*ca*nf*z3
       + 9620/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*ca*nf
       - 112/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*nf^2
       - 348*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf^2*nf
       + 112/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*ca*nf^2
       + 96*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*ca^2*nf*z3
       + 1520/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(2 + N)*S(R(1,-3),2 + N)*cf*ca*nf
       - 128*theta( - 3 + N)*den(2 + N)*S(R(1,-3),2 + N)*cf^2*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(1,-3),2 + N)*ca^2*nf
       - 80*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*cf^2*nf
       - 640/3*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*ca^2*nf
       + 192*theta( - 3 + N)*den(2 + N)*S(R(1,-2,1),2 + N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(2 + N)*S(R(1,-2,1),2 + N)*ca^2*nf
       + 4096/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*ca*nf
       + 112/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*nf^2
       - 244*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf^2*nf
       - 112/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*ca*nf^2
       - 188*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*ca^2*nf
       + 192*theta( - 3 + N)*den(2 + N)*S(R(1,1,-2),2 + N)*cf*ca*nf
       - 192*theta( - 3 + N)*den(2 + N)*S(R(1,1,-2),2 + N)*ca^2*nf
       + 416/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf*ca*nf
       - 16/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf*nf^2
       - 96*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf^2*nf
       + 16/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*ca*nf^2
       - 128/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*ca^2*nf
       + 128*theta( - 3 + N)*den(2 + N)*S(R(1,1,1,1),2 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(1,1,1,1),2 + N)*cf^2*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(1,1,1,1),2 + N)*ca^2*nf
       + 96*theta( - 3 + N)*den(2 + N)*S(R(1,1,2),2 + N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(2 + N)*S(R(1,1,2),2 + N)*cf^2*nf
       - 128*theta( - 3 + N)*den(2 + N)*S(R(1,1,2),2 + N)*ca^2*nf
       + 280/3*theta( - 3 + N)*den(2 + N)*S(R(1,2),2 + N)*cf*ca*nf
       + 80*theta( - 3 + N)*den(2 + N)*S(R(1,2),2 + N)*cf^2*nf
       - 16/3*theta( - 3 + N)*den(2 + N)*S(R(1,2),2 + N)*ca*nf^2
       - 168*theta( - 3 + N)*den(2 + N)*S(R(1,2),2 + N)*ca^2*nf
       - 96*theta( - 3 + N)*den(2 + N)*S(R(1,2,1),2 + N)*cf*ca*nf
       + 160*theta( - 3 + N)*den(2 + N)*S(R(1,2,1),2 + N)*cf^2*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(1,2,1),2 + N)*ca^2*nf
       - 32*theta( - 3 + N)*den(2 + N)*S(R(1,3),2 + N)*cf*ca*nf
       - 128*theta( - 3 + N)*den(2 + N)*S(R(1,3),2 + N)*cf^2*nf
       + 32*theta( - 3 + N)*den(2 + N)*S(R(1,3),2 + N)*ca^2*nf
       - 212*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf*ca*nf
       - 128/9*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf*nf^2
       + 272*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf^2*nf
       - 276*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*ca^2*nf
       - 96*theta( - 3 + N)*den(2 + N)*S(R(2,-2),2 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(2 + N)*S(R(2,-2),2 + N)*cf^2*nf
       + 32*theta( - 3 + N)*den(2 + N)*S(R(2,-2),2 + N)*ca^2*nf
       - 120*theta( - 3 + N)*den(2 + N)*S(R(2,1),2 + N)*cf*ca*nf
       + 176*theta( - 3 + N)*den(2 + N)*S(R(2,1),2 + N)*cf^2*nf
       + 16/3*theta( - 3 + N)*den(2 + N)*S(R(2,1),2 + N)*ca*nf^2
       - 184/3*theta( - 3 + N)*den(2 + N)*S(R(2,1),2 + N)*ca^2*nf
       - 192*theta( - 3 + N)*den(2 + N)*S(R(2,1,1),2 + N)*cf*ca*nf
       + 192*theta( - 3 + N)*den(2 + N)*S(R(2,1,1),2 + N)*cf^2*nf
       - 192*theta( - 3 + N)*den(2 + N)*S(R(2,2),2 + N)*cf^2*nf
       + 64*theta( - 3 + N)*den(2 + N)*S(R(2,2),2 + N)*ca^2*nf
       + 32/3*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*cf*nf^2
       - 80*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*cf^2*nf
       + 16/3*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*ca*nf^2
       - 592/3*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*ca^2*nf
       + 192*theta( - 3 + N)*den(2 + N)*S(R(3,1),2 + N)*cf*ca*nf
       - 224*theta( - 3 + N)*den(2 + N)*S(R(3,1),2 + N)*cf^2*nf
       - 96*theta( - 3 + N)*den(2 + N)*S(R(3,1),2 + N)*ca^2*nf
       - 48*theta( - 3 + N)*den(2 + N)*S(R(4),2 + N)*cf*ca*nf
       + 96*theta( - 3 + N)*den(2 + N)*S(R(4),2 + N)*cf^2*nf
       + 80*theta( - 3 + N)*den(2 + N)*S(R(4),2 + N)*ca^2*nf
       + 384*theta( - 3 + N)*den(N)*cf*ca*nf*z3
       - 1429/12*theta( - 3 + N)*den(N)*cf*ca*nf
       + 13025/54*theta( - 3 + N)*den(N)*cf*nf^2
       - 216*theta( - 3 + N)*den(N)*cf^2*nf*z3
       - 963/4*theta( - 3 + N)*den(N)*cf^2*nf
       - 226/27*theta( - 3 + N)*den(N)*ca*nf^2
       - 168*theta( - 3 + N)*den(N)*ca^2*nf*z3
       - 10454/27*theta( - 3 + N)*den(N)*ca^2*nf
       - 192*theta( - 3 + N)*den(N)^2*cf*ca*nf*z3
       + 2696/27*theta( - 3 + N)*den(N)^2*cf*ca*nf
       - 4856/27*theta( - 3 + N)*den(N)^2*cf*nf^2
       + 48*theta( - 3 + N)*den(N)^2*cf^2*nf*z3
       + 170*theta( - 3 + N)*den(N)^2*cf^2*nf
       - 160/27*theta( - 3 + N)*den(N)^2*ca*nf^2
       + 144*theta( - 3 + N)*den(N)^2*ca^2*nf*z3
       - 8696/27*theta( - 3 + N)*den(N)^2*ca^2*nf
       + 1267/9*theta( - 3 + N)*den(N)^3*cf*ca*nf
       + 1370/9*theta( - 3 + N)*den(N)^3*cf*nf^2
       - 89*theta( - 3 + N)*den(N)^3*cf^2*nf
       - 64/9*theta( - 3 + N)*den(N)^3*ca*nf^2
       - 2396/9*theta( - 3 + N)*den(N)^3*ca^2*nf
       + 44/3*theta( - 3 + N)*den(N)^4*cf*ca*nf
       - 152/3*theta( - 3 + N)*den(N)^4*cf*nf^2
       + 12*theta( - 3 + N)*den(N)^4*cf^2*nf
       - 16/3*theta( - 3 + N)*den(N)^4*ca*nf^2
       - 56/3*theta( - 3 + N)*den(N)^4*ca^2*nf
       + 32*theta( - 3 + N)*den(N)^5*cf*ca*nf
       + 32*theta( - 3 + N)*den(N)^5*cf*nf^2
       - 16*theta( - 3 + N)*den(N)^5*cf^2*nf
       - 64*theta( - 3 + N)*den(N)^5*ca^2*nf
       - 96*theta( - 3 + N)*den(N)^4*S(R(1),N)*cf*ca*nf
       + 56*theta( - 3 + N)*den(N)^4*S(R(1),N)*cf^2*nf
       + 12*theta( - 3 + N)*den(N)^4*S(R(1),N)*ca^2*nf
       - 48*theta( - 3 + N)*den(N)^3*S(R(-2),N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(N)^3*S(R(-2),N)*cf^2*nf
       - 72*theta( - 3 + N)*den(N)^3*S(R(-2),N)*ca^2*nf
       - 40*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf*ca*nf
       - 36*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf^2*nf
       - 44*theta( - 3 + N)*den(N)^3*S(R(1),N)*ca^2*nf
       + 64*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf*ca*nf
       - 80*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf^2*nf
       - 48*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*ca^2*nf
       - 64*theta( - 3 + N)*den(N)^3*S(R(2),N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(N)^3*S(R(2),N)*cf^2*nf
       + 64*theta( - 3 + N)*den(N)^2*S(R(-3),N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(N)^2*S(R(-3),N)*cf^2*nf
       - 24*theta( - 3 + N)*den(N)^2*S(R(-3),N)*ca^2*nf
       + 208*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf*ca*nf
       - 176*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf^2*nf
       - 72*theta( - 3 + N)*den(N)^2*S(R(-2),N)*ca^2*nf
       - 96*theta( - 3 + N)*den(N)^2*S(R(-2,1),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(-2,1),N)*ca^2*nf
       - 688/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*ca*nf
       - 80/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*nf^2
       + 68*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf^2*nf
       + 628/3*theta( - 3 + N)*den(N)^2*S(R(1),N)*ca^2*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*cf^2*nf
       - 48*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*ca^2*nf
       - 44/3*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf*ca*nf
       + 8/3*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf*nf^2
       + 36*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf^2*nf
       - 56*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf*ca*nf
       + 40*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf^2*nf
       + 16*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*ca^2*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*cf^2*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*ca^2*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf*ca*nf
       - 48*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf^2*nf
       - 16*theta( - 3 + N)*den(N)^2*S(R(2),N)*ca^2*nf
       - 64*theta( - 3 + N)*den(N)^2*S(R(2,1),N)*cf^2*nf
       + 88*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf^2*nf
       - 84*theta( - 3 + N)*den(N)^2*S(R(3),N)*ca^2*nf
       + 88*theta( - 3 + N)*den(N)*S(R(-4),N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(N)*S(R(-4),N)*cf^2*nf
       + 24*theta( - 3 + N)*den(N)*S(R(-4),N)*ca^2*nf
       + 40*theta( - 3 + N)*den(N)*S(R(-3),N)*cf*ca*nf
       - 128*theta( - 3 + N)*den(N)*S(R(-3),N)*cf^2*nf
       + 16/3*theta( - 3 + N)*den(N)*S(R(-3),N)*ca*nf^2
       - 256/3*theta( - 3 + N)*den(N)*S(R(-3),N)*ca^2*nf
       - 96*theta( - 3 + N)*den(N)*S(R(-3,1),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)*S(R(-3,1),N)*ca^2*nf
       - 296*theta( - 3 + N)*den(N)*S(R(-2),N)*cf*ca*nf
       + 320*theta( - 3 + N)*den(N)*S(R(-2),N)*cf^2*nf
       - 80/9*theta( - 3 + N)*den(N)*S(R(-2),N)*ca*nf^2
       + 704/9*theta( - 3 + N)*den(N)*S(R(-2),N)*ca^2*nf
       + 48*theta( - 3 + N)*den(N)*S(R(-2,-2),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)*S(R(-2,-2),N)*cf^2*nf
       - 48*theta( - 3 + N)*den(N)*S(R(-2,-2),N)*ca^2*nf
       + 112*theta( - 3 + N)*den(N)*S(R(-2,1),N)*cf*ca*nf
       + 64*theta( - 3 + N)*den(N)*S(R(-2,1),N)*ca^2*nf
       + 32*theta( - 3 + N)*den(N)*S(R(-2,1,1),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)*S(R(-2,1,1),N)*ca^2*nf
       - 32*theta( - 3 + N)*den(N)*S(R(-2,2),N)*cf*ca*nf
       - 48*theta( - 3 + N)*den(N)*S(R(1),N)*cf*ca*nf*z3
       + 3692/27*theta( - 3 + N)*den(N)*S(R(1),N)*cf*ca*nf
       - 224/27*theta( - 3 + N)*den(N)*S(R(1),N)*cf*nf^2
       + 32*theta( - 3 + N)*den(N)*S(R(1),N)*cf^2*nf
       + 152/27*theta( - 3 + N)*den(N)*S(R(1),N)*ca*nf^2
       + 48*theta( - 3 + N)*den(N)*S(R(1),N)*ca^2*nf*z3
       - 8288/27*theta( - 3 + N)*den(N)*S(R(1),N)*ca^2*nf
       - 96*theta( - 3 + N)*den(N)*S(R(1,-3),N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(N)*S(R(1,-3),N)*cf^2*nf
       - 32*theta( - 3 + N)*den(N)*S(R(1,-3),N)*ca^2*nf
       - 144*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf*ca*nf
       + 128*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf^2*nf
       + 112*theta( - 3 + N)*den(N)*S(R(1,-2),N)*ca^2*nf
       + 96*theta( - 3 + N)*den(N)*S(R(1,-2,1),N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(N)*S(R(1,-2,1),N)*ca^2*nf
       + 808/9*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*ca*nf
       + 80/9*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*nf^2
       - 108*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf^2*nf
       - 80/9*theta( - 3 + N)*den(N)*S(R(1,1),N)*ca*nf^2
       + 164/9*theta( - 3 + N)*den(N)*S(R(1,1),N)*ca^2*nf
       + 96*theta( - 3 + N)*den(N)*S(R(1,1,-2),N)*cf*ca*nf
       - 96*theta( - 3 + N)*den(N)*S(R(1,1,-2),N)*ca^2*nf
       + 56/3*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf*ca*nf
       - 8/3*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf*nf^2
       + 4*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf^2*nf
       + 8/3*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*ca*nf^2
       - 68/3*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*ca^2*nf
       + 64*theta( - 3 + N)*den(N)*S(R(1,1,1,1),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)*S(R(1,1,1,1),N)*cf^2*nf
       - 32*theta( - 3 + N)*den(N)*S(R(1,1,1,1),N)*ca^2*nf
       + 48*theta( - 3 + N)*den(N)*S(R(1,1,2),N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(N)*S(R(1,1,2),N)*cf^2*nf
       - 64*theta( - 3 + N)*den(N)*S(R(1,1,2),N)*ca^2*nf
       + 4*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf^2*nf
       - 8/3*theta( - 3 + N)*den(N)*S(R(1,2),N)*ca*nf^2
       + 92/3*theta( - 3 + N)*den(N)*S(R(1,2),N)*ca^2*nf
       - 48*theta( - 3 + N)*den(N)*S(R(1,2,1),N)*cf*ca*nf
       + 80*theta( - 3 + N)*den(N)*S(R(1,2,1),N)*cf^2*nf
       - 32*theta( - 3 + N)*den(N)*S(R(1,2,1),N)*ca^2*nf
       - 16*theta( - 3 + N)*den(N)*S(R(1,3),N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(N)*S(R(1,3),N)*cf^2*nf
       + 16*theta( - 3 + N)*den(N)*S(R(1,3),N)*ca^2*nf
       - 116/9*theta( - 3 + N)*den(N)*S(R(2),N)*cf*ca*nf
       - 40/9*theta( - 3 + N)*den(N)*S(R(2),N)*cf*nf^2
       + 160*theta( - 3 + N)*den(N)*S(R(2),N)*cf^2*nf
       + 92/3*theta( - 3 + N)*den(N)*S(R(2),N)*ca^2*nf
       - 48*theta( - 3 + N)*den(N)*S(R(2,-2),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)*S(R(2,-2),N)*cf^2*nf
       + 16*theta( - 3 + N)*den(N)*S(R(2,-2),N)*ca^2*nf
       - 20*theta( - 3 + N)*den(N)*S(R(2,1),N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(N)*S(R(2,1),N)*cf^2*nf
       + 8/3*theta( - 3 + N)*den(N)*S(R(2,1),N)*ca*nf^2
       - 44/3*theta( - 3 + N)*den(N)*S(R(2,1),N)*ca^2*nf
       - 96*theta( - 3 + N)*den(N)*S(R(2,1,1),N)*cf*ca*nf
       + 96*theta( - 3 + N)*den(N)*S(R(2,1,1),N)*cf^2*nf
       - 96*theta( - 3 + N)*den(N)*S(R(2,2),N)*cf^2*nf
       + 32*theta( - 3 + N)*den(N)*S(R(2,2),N)*ca^2*nf
       - 472/3*theta( - 3 + N)*den(N)*S(R(3),N)*cf*ca*nf
       + 16/3*theta( - 3 + N)*den(N)*S(R(3),N)*cf*nf^2
       + 32*theta( - 3 + N)*den(N)*S(R(3),N)*cf^2*nf
       + 8/3*theta( - 3 + N)*den(N)*S(R(3),N)*ca*nf^2
       + 304/3*theta( - 3 + N)*den(N)*S(R(3),N)*ca^2*nf
       + 96*theta( - 3 + N)*den(N)*S(R(3,1),N)*cf*ca*nf
       - 112*theta( - 3 + N)*den(N)*S(R(3,1),N)*cf^2*nf
       - 48*theta( - 3 + N)*den(N)*S(R(3,1),N)*ca^2*nf
       - 24*theta( - 3 + N)*den(N)*S(R(4),N)*cf*ca*nf
       + 48*theta( - 3 + N)*den(N)*S(R(4),N)*cf^2*nf
       + 40*theta( - 3 + N)*den(N)*S(R(4),N)*ca^2*nf
       - 104/3*delta( - 2 + N)*cf*ca*nf*z3
       + 139/9*delta( - 2 + N)*cf*ca*nf
       - 173/243*delta( - 2 + N)*cf*nf^2
       + 32/3*delta( - 2 + N)*cf^2*nf*z3
       - 2155/243*delta( - 2 + N)*cf^2*nf
       + 1058/243*delta( - 2 + N)*ca*nf^2
       + 24*delta( - 2 + N)*ca^2*nf*z3
       - 3589/162*delta( - 2 + N)*ca^2*nf
      ;


L   ggq2 =
       + 64*theta( - 3 + N)*den( - 1 + N)*cf*ca*nf*z3
       - 1934/9*theta( - 3 + N)*den( - 1 + N)*cf*ca*nf
       - 138305/81*theta( - 3 + N)*den( - 1 + N)*cf*ca^2
       + 16/9*theta( - 3 + N)*den( - 1 + N)*cf*nf^2
       - 163*theta( - 3 + N)*den( - 1 + N)*cf^2*ca
       - 64*theta( - 3 + N)*den( - 1 + N)*cf^2*nf*z3
       + 3250/9*theta( - 3 + N)*den( - 1 + N)*cf^2*nf
       + 94*theta( - 3 + N)*den( - 1 + N)*cf^3
       + 1208/27*theta( - 3 + N)*den( - 1 + N)^2*cf*ca*nf
       + 6320/27*theta( - 3 + N)*den( - 1 + N)^2*cf*ca^2
       - 1520/27*theta( - 3 + N)*den( - 1 + N)^2*cf^2*nf
       - 128*theta( - 3 + N)*den( - 1 + N)^3*S(R(1,1), - 1 + N)*cf*ca^2
       + 128*theta( - 3 + N)*den( - 1 + N)^3*S(R(2), - 1 + N)*cf*ca^2
       + 224*theta( - 3 + N)*den( - 1 + N)^2*S(R(-3), - 1 + N)*cf*ca^2
       + 32/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(-2), - 1 + N)*cf*ca*nf
       + 176/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(-2), - 1 + N)*cf*ca^2
       - 64/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(-2), - 1 + N)*cf^2*nf
       - 128*theta( - 3 + N)*den( - 1 + N)^2*S(R(-2,1), - 1 + N)*cf*ca^2
       - 320/9*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*cf*ca*nf
       + 128*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*cf^2*ca
       - 64*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,-2), - 1 + N)*cf*ca^2
       + 80/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,1), - 1 + N)*cf*ca*nf
       - 544/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,1), - 1 + N)*cf*ca^2
       + 144*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,1), - 1 + N)*cf^2*ca
       + 64*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,1,1), - 1 + N)*cf*ca^2
       - 64*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,1,1), - 1 + N)*cf^2*ca
       - 96*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,2), - 1 + N)*cf*ca^2
       - 32*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,2), - 1 + N)*cf^2*ca
       - 32/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(2), - 1 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(2), - 1 + N)*cf*ca^2
       - 96*theta( - 3 + N)*den( - 1 + N)^2*S(R(2), - 1 + N)*cf^2*ca
       - 96*theta( - 3 + N)*den( - 1 + N)^2*S(R(2,1), - 1 + N)*cf*ca^2
       - 32*theta( - 3 + N)*den( - 1 + N)^2*S(R(2,1), - 1 + N)*cf^2*ca
       + 96*theta( - 3 + N)*den( - 1 + N)^2*S(R(3), - 1 + N)*cf*ca^2
       + 64*theta( - 3 + N)*den( - 1 + N)^2*S(R(3), - 1 + N)*cf^2*ca
       + 176*theta( - 3 + N)*den( - 1 + N)*S(R(-4), - 1 + N)*cf*ca^2
       + 48*theta( - 3 + N)*den( - 1 + N)*S(R(-4), - 1 + N)*cf^2*ca
       + 32*theta( - 3 + N)*den( - 1 + N)*S(R(-4), - 1 + N)*cf^3
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*cf*ca*nf
       - 784/3*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*cf*ca^2
       + 48*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*cf^2*ca
       - 192*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*cf^3
       - 256*theta( - 3 + N)*den( - 1 + N)*S(R(-3,1), - 1 + N)*cf*ca^2
       - 656/9*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*cf*ca*nf
       - 872/3*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*cf*ca^2
       - 32*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*cf^2*ca
       + 496/9*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*cf^2*nf
       + 112*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*cf^3
       - 32*theta( - 3 + N)*den( - 1 + N)*S(R(-2,-2), - 1 + N)*cf*ca^2
       + 160*theta( - 3 + N)*den( - 1 + N)*S(R(-2,-2), - 1 + N)*cf^2*ca
       - 192*theta( - 3 + N)*den( - 1 + N)*S(R(-2,-2), - 1 + N)*cf^3
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1), - 1 + N)*cf*ca*nf
       + 48*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1), - 1 + N)*cf*ca^2
       + 96*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1), - 1 + N)*cf^2*ca
       + 64*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1,1), - 1 + N)*cf*ca^2
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1,1), - 1 + N)*cf^2*ca
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(-2,2), - 1 + N)*cf*ca^2
       - 2008/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*ca*nf
       - 96*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*ca^2*z3
       + 13516/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*ca^2
       + 80/9*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*nf^2
       + 288*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf^2*ca*z3
       - 13006/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf^2*ca
       + 56/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf^2*nf
       - 192*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf^3*z3
       + 94*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf^3
       - 384*theta( - 3 + N)*den( - 1 + N)*S(R(1,-3), - 1 + N)*cf*ca^2
       + 128*theta( - 3 + N)*den( - 1 + N)*S(R(1,-3), - 1 + N)*cf^2*ca
       - 128*theta( - 3 + N)*den( - 1 + N)*S(R(1,-3), - 1 + N)*cf^3
       + 128/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*cf*ca*nf
       - 304*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*cf*ca^2
       - 96*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*cf^2*ca
       + 192*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*cf^3
       + 320*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2,1), - 1 + N)*cf*ca^2
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2,1), - 1 + N)*cf^2*ca
       + 104*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*ca*nf
       - 260/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*ca^2
       - 16/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*nf^2
       - 64/9*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf^2*ca
       - 664/9*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf^2*nf
       + 92*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf^3
       + 192*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,-2), - 1 + N)*cf*ca^2
       - 192*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,-2), - 1 + N)*cf^2*ca
       - 80/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf*ca*nf
       + 688/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf*ca^2
       - 1120/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf^2*ca
       + 80/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf^2*nf
       + 144*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf^3
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1,1), - 1 + N)*cf*ca^2
       + 128*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1,1), - 1 + N)*cf^2*ca
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1,1), - 1 + N)*cf^3
       + 128*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,2), - 1 + N)*cf*ca^2
       - 96*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,2), - 1 + N)*cf^2*ca
       - 32*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,2), - 1 + N)*cf^3
       + 80/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf*ca*nf
       - 440/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf*ca^2
       + 8/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf^2*ca
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf^2*nf
       + 96*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf^3
       + 128*theta( - 3 + N)*den( - 1 + N)*S(R(1,2,1), - 1 + N)*cf*ca^2
       - 32*theta( - 3 + N)*den( - 1 + N)*S(R(1,2,1), - 1 + N)*cf^2*ca
       - 96*theta( - 3 + N)*den( - 1 + N)*S(R(1,2,1), - 1 + N)*cf^3
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(1,3), - 1 + N)*cf*ca^2
       - 160*theta( - 3 + N)*den( - 1 + N)*S(R(1,3), - 1 + N)*cf^2*ca
       + 96*theta( - 3 + N)*den( - 1 + N)*S(R(1,3), - 1 + N)*cf^3
       - 320/9*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf*ca*nf
       - 284/3*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf*ca^2
       + 764/3*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf^2*ca
       - 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf^2*nf
       - 56*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf^3
       - 160*theta( - 3 + N)*den( - 1 + N)*S(R(2,-2), - 1 + N)*cf*ca^2
       + 96*theta( - 3 + N)*den( - 1 + N)*S(R(2,-2), - 1 + N)*cf^2*ca
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(2,-2), - 1 + N)*cf^3
       + 80/3*theta( - 3 + N)*den( - 1 + N)*S(R(2,1), - 1 + N)*cf*ca*nf
       - 440/3*theta( - 3 + N)*den( - 1 + N)*S(R(2,1), - 1 + N)*cf*ca^2
       + 328/3*theta( - 3 + N)*den( - 1 + N)*S(R(2,1), - 1 + N)*cf^2*ca
       + 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(2,1), - 1 + N)*cf^2*nf
       + 128*theta( - 3 + N)*den( - 1 + N)*S(R(2,1,1), - 1 + N)*cf*ca^2
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(2,1,1), - 1 + N)*cf^2*ca
       - 64*theta( - 3 + N)*den( - 1 + N)*S(R(2,1,1), - 1 + N)*cf^3
       - 192*theta( - 3 + N)*den( - 1 + N)*S(R(2,2), - 1 + N)*cf*ca^2
       + 64*theta( - 3 + N)*den( - 1 + N)*S(R(2,2), - 1 + N)*cf^3
       + 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf*ca*nf
       - 176*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf*ca^2
       - 128*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf^2*ca
       - 16*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf^2*nf
       + 48*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf^3
       - 128*theta( - 3 + N)*den( - 1 + N)*S(R(3,1), - 1 + N)*cf*ca^2
       + 48*theta( - 3 + N)*den( - 1 + N)*S(R(4), - 1 + N)*cf*ca^2
       + 80*theta( - 3 + N)*den( - 1 + N)*S(R(4), - 1 + N)*cf^2*ca
       + 32*theta( - 3 + N)*den(1 + N)*cf*ca*nf*z3
       - 148/3*theta( - 3 + N)*den(1 + N)*cf*ca*nf
       - 120*theta( - 3 + N)*den(1 + N)*cf*ca^2*z3
       - 1789/3*theta( - 3 + N)*den(1 + N)*cf*ca^2
       - 32/9*theta( - 3 + N)*den(1 + N)*cf*nf^2
       + 360*theta( - 3 + N)*den(1 + N)*cf^2*ca*z3
       - 3589/36*theta( - 3 + N)*den(1 + N)*cf^2*ca
       - 32*theta( - 3 + N)*den(1 + N)*cf^2*nf*z3
       - 3301/54*theta( - 3 + N)*den(1 + N)*cf^2*nf
       - 240*theta( - 3 + N)*den(1 + N)*cf^3*z3
       + 363/4*theta( - 3 + N)*den(1 + N)*cf^3
       + 1688/27*theta( - 3 + N)*den(1 + N)^2*cf*ca*nf
       + 3202/27*theta( - 3 + N)*den(1 + N)^2*cf*ca^2
       + 7061/27*theta( - 3 + N)*den(1 + N)^2*cf^2*ca
       - 3134/27*theta( - 3 + N)*den(1 + N)^2*cf^2*nf
       - 245*theta( - 3 + N)*den(1 + N)^2*cf^3
       + 668/9*theta( - 3 + N)*den(1 + N)^3*cf*ca*nf
       - 392/9*theta( - 3 + N)*den(1 + N)^3*cf*ca^2
       - 4631/9*theta( - 3 + N)*den(1 + N)^3*cf^2*ca
       - 1306/9*theta( - 3 + N)*den(1 + N)^3*cf^2*nf
       + 555*theta( - 3 + N)*den(1 + N)^3*cf^3
       - 224*theta( - 3 + N)*den(1 + N)^4*cf*ca^2
       + 608/3*theta( - 3 + N)*den(1 + N)^4*cf^2*ca
       - 32/3*theta( - 3 + N)*den(1 + N)^4*cf^2*nf
       - 132*theta( - 3 + N)*den(1 + N)^4*cf^3
       - 120*theta( - 3 + N)*den(1 + N)^5*cf*ca^2
       - 64*theta( - 3 + N)*den(1 + N)^5*cf^2*ca
       + 32*theta( - 3 + N)*den(1 + N)^5*cf^2*nf
       + 32*theta( - 3 + N)*den(1 + N)^5*cf^3
       + 100*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf*ca^2
       + 32*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf^2*ca
       + 24*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf^3
       + 8*theta( - 3 + N)*den(1 + N)^3*S(R(-2),1 + N)*cf*ca^2
       + 16*theta( - 3 + N)*den(1 + N)^3*S(R(-2),1 + N)*cf^2*ca
       - 32/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf*ca*nf
       + 980/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf*ca^2
       - 392/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf^2*ca
       - 16/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf^2*nf
       - 44*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf^3
       - 96*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf*ca^2
       + 8*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf^2*ca
       + 24*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf^3
       + 128*theta( - 3 + N)*den(1 + N)^3*S(R(2),1 + N)*cf*ca^2
       + 32*theta( - 3 + N)*den(1 + N)^3*S(R(2),1 + N)*cf^2*ca
       + 152*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*cf*ca^2
       + 32*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*cf^3
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf*ca*nf
       - 16*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf*ca^2
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf^2*ca
       - 32*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf^2*nf
       - 80*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf^3
       - 64*theta( - 3 + N)*den(1 + N)^2*S(R(-2,1),1 + N)*cf*ca^2
       - 236/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*ca*nf
       + 68/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*ca^2
       - 944/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf^2*ca
       + 128/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf^2*nf
       + 88*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf^3
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(1,-2),1 + N)*cf*ca^2
       + 8*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf*ca*nf
       - 132*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf*ca^2
       + 448/3*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf^2*ca
       + 8/3*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf^2*nf
       - 52*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf^3
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf*ca^2
       - 8*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf^2*ca
       - 8*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf^3
       - 48*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*cf*ca^2
       - 16*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*cf^2*ca
       - 16/3*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf*ca*nf
       + 592/3*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf*ca^2
       - 224/3*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf^2*ca
       - 16/3*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf^2*nf
       - 48*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf^3
       - 48*theta( - 3 + N)*den(1 + N)^2*S(R(2,1),1 + N)*cf*ca^2
       - 32*theta( - 3 + N)*den(1 + N)^2*S(R(2,1),1 + N)*cf^2*ca
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(2,1),1 + N)*cf^3
       + 84*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf*ca^2
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf^2*ca
       - 8*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf^3
       + 88*theta( - 3 + N)*den(1 + N)*S(R(-4),1 + N)*cf*ca^2
       + 24*theta( - 3 + N)*den(1 + N)*S(R(-4),1 + N)*cf^2*ca
       + 16*theta( - 3 + N)*den(1 + N)*S(R(-4),1 + N)*cf^3
       - 32/3*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf*ca*nf
       - 232/3*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf*ca^2
       + 144*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf^2*ca
       - 224*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf^3
       - 128*theta( - 3 + N)*den(1 + N)*S(R(-3,1),1 + N)*cf*ca^2
       + 8*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf*ca*nf
       + 212/3*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf*ca^2
       - 700*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf^2*ca
       - 48*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf^2*nf
       + 608*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf^3
       - 16*theta( - 3 + N)*den(1 + N)*S(R(-2,-2),1 + N)*cf*ca^2
       + 80*theta( - 3 + N)*den(1 + N)*S(R(-2,-2),1 + N)*cf^2*ca
       - 96*theta( - 3 + N)*den(1 + N)*S(R(-2,-2),1 + N)*cf^3
       + 32/3*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*cf*ca*nf
       + 304/3*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*cf*ca^2
       + 64*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*cf^2*ca
       + 32*theta( - 3 + N)*den(1 + N)*S(R(-2,1,1),1 + N)*cf*ca^2
       - 32*theta( - 3 + N)*den(1 + N)*S(R(-2,1,1),1 + N)*cf^2*ca
       - 32*theta( - 3 + N)*den(1 + N)*S(R(-2,2),1 + N)*cf*ca^2
       - 1720/27*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*ca*nf
       - 48*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*ca^2*z3
       + 766/27*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*ca^2
       + 64/9*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*nf^2
       + 144*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^2*ca*z3
       - 10936/27*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^2*ca
       + 1360/27*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^2*nf
       - 96*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^3*z3
       + 258*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^3
       - 192*theta( - 3 + N)*den(1 + N)*S(R(1,-3),1 + N)*cf*ca^2
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,-3),1 + N)*cf^2*ca
       - 64*theta( - 3 + N)*den(1 + N)*S(R(1,-3),1 + N)*cf^3
       + 64/3*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf*ca*nf
       + 176/3*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf*ca^2
       - 304*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf^2*ca
       + 320*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf^3
       + 160*theta( - 3 + N)*den(1 + N)*S(R(1,-2,1),1 + N)*cf*ca^2
       - 32*theta( - 3 + N)*den(1 + N)*S(R(1,-2,1),1 + N)*cf^2*ca
       + 488/9*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*ca*nf
       - 914/9*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*ca^2
       - 8/3*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*nf^2
       + 1448/9*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf^2*ca
       - 512/9*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf^2*nf
       - 54*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf^3
       + 96*theta( - 3 + N)*den(1 + N)*S(R(1,1,-2),1 + N)*cf*ca^2
       - 96*theta( - 3 + N)*den(1 + N)*S(R(1,1,-2),1 + N)*cf^2*ca
       - 40/3*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf*ca*nf
       + 340/3*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf*ca^2
       - 712/3*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf^2*ca
       + 40/3*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf^2*nf
       + 124*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf^3
       - 32*theta( - 3 + N)*den(1 + N)*S(R(1,1,1,1),1 + N)*cf*ca^2
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,1,1,1),1 + N)*cf^2*ca
       - 32*theta( - 3 + N)*den(1 + N)*S(R(1,1,1,1),1 + N)*cf^3
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,1,2),1 + N)*cf*ca^2
       - 48*theta( - 3 + N)*den(1 + N)*S(R(1,1,2),1 + N)*cf^2*ca
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1,1,2),1 + N)*cf^3
       + 40/3*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf*ca*nf
       - 412/3*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf*ca^2
       + 148/3*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf^2*ca
       + 32/3*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf^2*nf
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf^3
       + 64*theta( - 3 + N)*den(1 + N)*S(R(1,2,1),1 + N)*cf*ca^2
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1,2,1),1 + N)*cf^2*ca
       - 48*theta( - 3 + N)*den(1 + N)*S(R(1,2,1),1 + N)*cf^3
       - 32*theta( - 3 + N)*den(1 + N)*S(R(1,3),1 + N)*cf*ca^2
       - 80*theta( - 3 + N)*den(1 + N)*S(R(1,3),1 + N)*cf^2*ca
       + 48*theta( - 3 + N)*den(1 + N)*S(R(1,3),1 + N)*cf^3
       - 280/9*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*ca*nf
       + 400/9*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*ca^2
       + 80/3*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf^2*ca
       - 8/3*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf^2*nf
       - 24*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf^3
       - 80*theta( - 3 + N)*den(1 + N)*S(R(2,-2),1 + N)*cf*ca^2
       + 48*theta( - 3 + N)*den(1 + N)*S(R(2,-2),1 + N)*cf^2*ca
       - 32*theta( - 3 + N)*den(1 + N)*S(R(2,-2),1 + N)*cf^3
       + 40/3*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*cf*ca*nf
       - 412/3*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*cf*ca^2
       + 332/3*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*cf^2*ca
       + 16/3*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*cf^2*nf
       + 8*theta( - 3 + N)*den(1 + N)*S(R(2,1),1 + N)*cf^3
       + 64*theta( - 3 + N)*den(1 + N)*S(R(2,1,1),1 + N)*cf*ca^2
       - 32*theta( - 3 + N)*den(1 + N)*S(R(2,1,1),1 + N)*cf^2*ca
       - 32*theta( - 3 + N)*den(1 + N)*S(R(2,1,1),1 + N)*cf^3
       - 96*theta( - 3 + N)*den(1 + N)*S(R(2,2),1 + N)*cf*ca^2
       + 32*theta( - 3 + N)*den(1 + N)*S(R(2,2),1 + N)*cf^3
       + 16/3*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf*ca*nf
       + 332/3*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf*ca^2
       - 232*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf^2*ca
       - 8*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf^2*nf
       + 108*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf^3
       - 64*theta( - 3 + N)*den(1 + N)*S(R(3,1),1 + N)*cf*ca^2
       + 24*theta( - 3 + N)*den(1 + N)*S(R(4),1 + N)*cf*ca^2
       + 40*theta( - 3 + N)*den(1 + N)*S(R(4),1 + N)*cf^2*ca
       + 1048/9*theta( - 3 + N)*den(2 + N)*cf*ca*nf
       + 33680/81*theta( - 3 + N)*den(2 + N)*cf*ca^2
       - 56*theta( - 3 + N)*den(2 + N)*cf^2*ca
       + 224/9*theta( - 3 + N)*den(2 + N)*cf^2*nf
       + 2000/27*theta( - 3 + N)*den(2 + N)^2*cf*ca*nf
       + 9344/27*theta( - 3 + N)*den(2 + N)^2*cf*ca^2
       + 2680/27*theta( - 3 + N)*den(2 + N)^2*cf^2*ca
       - 352/27*theta( - 3 + N)*den(2 + N)^2*cf^2*nf
       + 32/9*theta( - 3 + N)*den(2 + N)^3*cf*ca*nf
       + 616/9*theta( - 3 + N)*den(2 + N)^3*cf*ca^2
       + 448/9*theta( - 3 + N)*den(2 + N)^3*cf^2*ca
       - 64/9*theta( - 3 + N)*den(2 + N)^3*cf^2*nf
       - 224/3*theta( - 3 + N)*den(2 + N)^4*cf*ca^2
       - 128/3*theta( - 3 + N)*den(2 + N)^4*cf^2*ca
       + 96*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*cf*ca^2
       + 64*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*cf*ca^2
       - 32/3*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*ca*nf
       + 1664/9*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*ca^2
       - 96*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf^2*ca
       - 128/3*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf*ca^2
       + 32*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf^2*ca
       + 224/3*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*cf*ca^2
       + 128/3*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*cf*ca^2
       + 32/9*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*cf*ca*nf
       + 160*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*cf*ca^2
       - 64/9*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*cf^2*nf
       - 64/3*theta( - 3 + N)*den(2 + N)*S(R(-2,1),2 + N)*cf*ca^2
       - 400/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*ca*nf
       + 1744/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*ca^2
       - 1384/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf^2*ca
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*cf*ca^2
       - 32/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*ca*nf
       - 592/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*ca^2
       + 416/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf^2*ca
       - 32/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf*ca^2
       + 32/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf^2*ca
       + 1160/9*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf*ca^2
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*cf*ca^2
       - 64*theta( - 3 + N)*den(N)*cf*ca*nf*z3
       + 2320/9*theta( - 3 + N)*den(N)*cf*ca*nf
       + 192*theta( - 3 + N)*den(N)*cf*ca^2*z3
       + 12710/9*theta( - 3 + N)*den(N)*cf*ca^2
       - 16/9*theta( - 3 + N)*den(N)*cf*nf^2
       - 576*theta( - 3 + N)*den(N)*cf^2*ca*z3
       + 11629/18*theta( - 3 + N)*den(N)*cf^2*ca
       + 64*theta( - 3 + N)*den(N)*cf^2*nf*z3
       - 8707/27*theta( - 3 + N)*den(N)*cf^2*nf
       + 384*theta( - 3 + N)*den(N)*cf^3*z3
       - 475/2*theta( - 3 + N)*den(N)*cf^3
       + 6656/27*theta( - 3 + N)*den(N)^2*cf*ca*nf
       - 96*theta( - 3 + N)*den(N)^2*cf*ca^2*z3
       + 38224/27*theta( - 3 + N)*den(N)^2*cf*ca^2
       + 288*theta( - 3 + N)*den(N)^2*cf^2*ca*z3
       - 9925/27*theta( - 3 + N)*den(N)^2*cf^2*ca
       - 5174/27*theta( - 3 + N)*den(N)^2*cf^2*nf
       - 192*theta( - 3 + N)*den(N)^2*cf^3*z3
       + 91*theta( - 3 + N)*den(N)^2*cf^3
       + 1088/9*theta( - 3 + N)*den(N)^3*cf*ca*nf
       + 5080/9*theta( - 3 + N)*den(N)^3*cf*ca^2
       + 2590/9*theta( - 3 + N)*den(N)^3*cf^2*ca
       - 2236/9*theta( - 3 + N)*den(N)^3*cf^2*nf
       + 70*theta( - 3 + N)*den(N)^3*cf^3
       + 32/3*theta( - 3 + N)*den(N)^4*cf*ca*nf
       - 272/3*theta( - 3 + N)*den(N)^4*cf*ca^2
       - 760/3*theta( - 3 + N)*den(N)^4*cf^2*ca
       + 16/3*theta( - 3 + N)*den(N)^4*cf^2*nf
       + 16*theta( - 3 + N)*den(N)^4*cf^3
       + 128*theta( - 3 + N)*den(N)^5*cf*ca^2
       - 64*theta( - 3 + N)*den(N)^5*cf^2*ca
       - 64*theta( - 3 + N)*den(N)^5*cf^2*nf
       + 32*theta( - 3 + N)*den(N)^5*cf^3
       + 104*theta( - 3 + N)*den(N)^4*S(R(1),N)*cf*ca^2
       - 48*theta( - 3 + N)*den(N)^4*S(R(1),N)*cf^3
       + 144*theta( - 3 + N)*den(N)^3*S(R(-2),N)*cf*ca^2
       + 224*theta( - 3 + N)*den(N)^3*S(R(-2),N)*cf^2*ca
       - 256*theta( - 3 + N)*den(N)^3*S(R(-2),N)*cf^3
       - 64/3*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf*ca*nf
       + 2140/3*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf*ca^2
       - 752/3*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf^2*ca
       + 32/3*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf^2*nf
       + 64*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf^3
       + 64*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf*ca^2
       + 112*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf^2*ca
       - 48*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf^3
       + 64*theta( - 3 + N)*den(N)^3*S(R(2),N)*cf*ca^2
       - 64*theta( - 3 + N)*den(N)^3*S(R(2),N)*cf^2*ca
       - 176*theta( - 3 + N)*den(N)^2*S(R(-3),N)*cf*ca^2
       + 96*theta( - 3 + N)*den(N)^2*S(R(-3),N)*cf^2*ca
       - 128*theta( - 3 + N)*den(N)^2*S(R(-3),N)*cf^3
       + 160/3*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf*ca*nf
       + 1160/3*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf*ca^2
       - 576*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf^2*ca
       - 64*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf^2*nf
       + 480*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf^3
       + 256*theta( - 3 + N)*den(N)^2*S(R(-2,1),N)*cf*ca^2
       + 256/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*ca*nf
       - 568/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*ca^2
       - 2456/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf^2*ca
       + 128/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf^2*nf
       + 164*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf^3
       + 352*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*cf*ca^2
       - 192*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*cf^2*ca
       + 128*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*cf^3
       - 48*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf*ca*nf
       + 104*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf*ca^2
       - 176/3*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf^2*ca
       - 16/3*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf^2*nf
       - 16*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf^3
       - 128*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf*ca^2
       + 112*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf^2*ca
       + 16*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf^3
       + 96*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*cf*ca^2
       + 32*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*cf^2*ca
       + 32/3*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf*ca*nf
       + 784/3*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf*ca^2
       + 16/3*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf^2*ca
       + 32/3*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf^2*nf
       + 96*theta( - 3 + N)*den(N)^2*S(R(2,1),N)*cf*ca^2
       + 64*theta( - 3 + N)*den(N)^2*S(R(2,1),N)*cf^2*ca
       - 32*theta( - 3 + N)*den(N)^2*S(R(2,1),N)*cf^3
       + 104*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf*ca^2
       - 224*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf^2*ca
       + 80*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf^3
       - 176*theta( - 3 + N)*den(N)*S(R(-4),N)*cf*ca^2
       - 48*theta( - 3 + N)*den(N)*S(R(-4),N)*cf^2*ca
       - 32*theta( - 3 + N)*den(N)*S(R(-4),N)*cf^3
       + 64/3*theta( - 3 + N)*den(N)*S(R(-3),N)*cf*ca*nf
       + 1160/3*theta( - 3 + N)*den(N)*S(R(-3),N)*cf*ca^2
       - 240*theta( - 3 + N)*den(N)*S(R(-3),N)*cf^2*ca
       + 384*theta( - 3 + N)*den(N)*S(R(-3),N)*cf^3
       + 256*theta( - 3 + N)*den(N)*S(R(-3,1),N)*cf*ca^2
       + 80/3*theta( - 3 + N)*den(N)*S(R(-2),N)*cf*ca*nf
       + 920/3*theta( - 3 + N)*den(N)*S(R(-2),N)*cf*ca^2
       + 608*theta( - 3 + N)*den(N)*S(R(-2),N)*cf^2*ca
       - 688*theta( - 3 + N)*den(N)*S(R(-2),N)*cf^3
       + 32*theta( - 3 + N)*den(N)*S(R(-2,-2),N)*cf*ca^2
       - 160*theta( - 3 + N)*den(N)*S(R(-2,-2),N)*cf^2*ca
       + 192*theta( - 3 + N)*den(N)*S(R(-2,-2),N)*cf^3
       - 64/3*theta( - 3 + N)*den(N)*S(R(-2,1),N)*cf*ca*nf
       - 608/3*theta( - 3 + N)*den(N)*S(R(-2,1),N)*cf*ca^2
       - 96*theta( - 3 + N)*den(N)*S(R(-2,1),N)*cf^2*ca
       - 64*theta( - 3 + N)*den(N)*S(R(-2,1,1),N)*cf*ca^2
       + 64*theta( - 3 + N)*den(N)*S(R(-2,1,1),N)*cf^2*ca
       + 64*theta( - 3 + N)*den(N)*S(R(-2,2),N)*cf*ca^2
       + 1184/27*theta( - 3 + N)*den(N)*S(R(1),N)*cf*ca*nf
       + 96*theta( - 3 + N)*den(N)*S(R(1),N)*cf*ca^2*z3
       - 6260/27*theta( - 3 + N)*den(N)*S(R(1),N)*cf*ca^2
       - 80/9*theta( - 3 + N)*den(N)*S(R(1),N)*cf*nf^2
       - 288*theta( - 3 + N)*den(N)*S(R(1),N)*cf^2*ca*z3
       + 16148/27*theta( - 3 + N)*den(N)*S(R(1),N)*cf^2*ca
       - 8/27*theta( - 3 + N)*den(N)*S(R(1),N)*cf^2*nf
       + 192*theta( - 3 + N)*den(N)*S(R(1),N)*cf^3*z3
       - 316*theta( - 3 + N)*den(N)*S(R(1),N)*cf^3
       + 384*theta( - 3 + N)*den(N)*S(R(1,-3),N)*cf*ca^2
       - 128*theta( - 3 + N)*den(N)*S(R(1,-3),N)*cf^2*ca
       + 128*theta( - 3 + N)*den(N)*S(R(1,-3),N)*cf^3
       - 128/3*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf*ca*nf
       + 176/3*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf*ca^2
       + 480*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf^2*ca
       - 448*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf^3
       - 320*theta( - 3 + N)*den(N)*S(R(1,-2,1),N)*cf*ca^2
       + 64*theta( - 3 + N)*den(N)*S(R(1,-2,1),N)*cf^2*ca
       - 880/9*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*ca*nf
       + 100/9*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*ca^2
       + 16/3*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*nf^2
       + 884/9*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf^2*ca
       + 712/9*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf^2*nf
       - 96*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf^3
       - 192*theta( - 3 + N)*den(N)*S(R(1,1,-2),N)*cf*ca^2
       + 192*theta( - 3 + N)*den(N)*S(R(1,1,-2),N)*cf^2*ca
       + 80/3*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf*ca*nf
       - 632/3*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf*ca^2
       + 1112/3*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf^2*ca
       - 80/3*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf^2*nf
       - 160*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf^3
       + 64*theta( - 3 + N)*den(N)*S(R(1,1,1,1),N)*cf*ca^2
       - 128*theta( - 3 + N)*den(N)*S(R(1,1,1,1),N)*cf^2*ca
       + 64*theta( - 3 + N)*den(N)*S(R(1,1,1,1),N)*cf^3
       - 128*theta( - 3 + N)*den(N)*S(R(1,1,2),N)*cf*ca^2
       + 96*theta( - 3 + N)*den(N)*S(R(1,1,2),N)*cf^2*ca
       + 32*theta( - 3 + N)*den(N)*S(R(1,1,2),N)*cf^3
       - 80/3*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf*ca*nf
       + 440/3*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf*ca^2
       - 8/3*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf^2*ca
       - 64/3*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf^2*nf
       - 96*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf^3
       - 128*theta( - 3 + N)*den(N)*S(R(1,2,1),N)*cf*ca^2
       + 32*theta( - 3 + N)*den(N)*S(R(1,2,1),N)*cf^2*ca
       + 96*theta( - 3 + N)*den(N)*S(R(1,2,1),N)*cf^3
       + 64*theta( - 3 + N)*den(N)*S(R(1,3),N)*cf*ca^2
       + 160*theta( - 3 + N)*den(N)*S(R(1,3),N)*cf^2*ca
       - 96*theta( - 3 + N)*den(N)*S(R(1,3),N)*cf^3
       + 320/9*theta( - 3 + N)*den(N)*S(R(2),N)*cf*ca*nf
       + 1096/9*theta( - 3 + N)*den(N)*S(R(2),N)*cf*ca^2
       - 224*theta( - 3 + N)*den(N)*S(R(2),N)*cf^2*ca
       - 16*theta( - 3 + N)*den(N)*S(R(2),N)*cf^3
       + 160*theta( - 3 + N)*den(N)*S(R(2,-2),N)*cf*ca^2
       - 96*theta( - 3 + N)*den(N)*S(R(2,-2),N)*cf^2*ca
       + 64*theta( - 3 + N)*den(N)*S(R(2,-2),N)*cf^3
       - 80/3*theta( - 3 + N)*den(N)*S(R(2,1),N)*cf*ca*nf
       + 440/3*theta( - 3 + N)*den(N)*S(R(2,1),N)*cf*ca^2
       - 424/3*theta( - 3 + N)*den(N)*S(R(2,1),N)*cf^2*ca
       - 32/3*theta( - 3 + N)*den(N)*S(R(2,1),N)*cf^2*nf
       + 32*theta( - 3 + N)*den(N)*S(R(2,1),N)*cf^3
       - 128*theta( - 3 + N)*den(N)*S(R(2,1,1),N)*cf*ca^2
       + 64*theta( - 3 + N)*den(N)*S(R(2,1,1),N)*cf^2*ca
       + 64*theta( - 3 + N)*den(N)*S(R(2,1,1),N)*cf^3
       + 192*theta( - 3 + N)*den(N)*S(R(2,2),N)*cf*ca^2
       - 64*theta( - 3 + N)*den(N)*S(R(2,2),N)*cf^3
       - 32/3*theta( - 3 + N)*den(N)*S(R(3),N)*cf*ca*nf
       + 140/3*theta( - 3 + N)*den(N)*S(R(3),N)*cf*ca^2
       + 384*theta( - 3 + N)*den(N)*S(R(3),N)*cf^2*ca
       + 16*theta( - 3 + N)*den(N)*S(R(3),N)*cf^2*nf
       - 192*theta( - 3 + N)*den(N)*S(R(3),N)*cf^3
       + 128*theta( - 3 + N)*den(N)*S(R(3,1),N)*cf*ca^2
       - 48*theta( - 3 + N)*den(N)*S(R(4),N)*cf*ca^2
       - 80*theta( - 3 + N)*den(N)*S(R(4),N)*cf^2*ca
       + 128/3*delta( - 2 + N)*cf*ca*nf*z3
       + 22/9*delta( - 2 + N)*cf*ca*nf
       - 64/3*delta( - 2 + N)*cf*ca^2*z3
       - 20920/243*delta( - 2 + N)*cf*ca^2
       + 284/81*delta( - 2 + N)*cf*nf^2
       + 64*delta( - 2 + N)*cf^2*ca*z3
       + 8528/243*delta( - 2 + N)*cf^2*ca
       - 128/3*delta( - 2 + N)*cf^2*nf*z3
       + 7094/243*delta( - 2 + N)*cf^2*nf
       - 128/3*delta( - 2 + N)*cf^3*z3
       + 560/243*delta( - 2 + N)*cf^3
      ;

L   gqqps2 =
       - 64*theta( - 3 + N)*den( - 1 + N)*cf*ca*nf*z3
       + 27044/81*theta( - 3 + N)*den( - 1 + N)*cf*ca*nf
       - 64/27*theta( - 3 + N)*den( - 1 + N)*cf*nf^2
       + 64*theta( - 3 + N)*den( - 1 + N)*cf^2*nf*z3
       - 220/3*theta( - 3 + N)*den( - 1 + N)*cf^2*nf
       - 896/27*theta( - 3 + N)*den( - 1 + N)^2*cf*ca*nf
       + 320/9*theta( - 3 + N)*den( - 1 + N)^2*S(R(1), - 1 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den( - 1 + N)^2*S(R(1,1), - 1 + N)*cf*ca*nf
       + 128/3*theta( - 3 + N)*den( - 1 + N)*S(R(-3), - 1 + N)*cf*ca*nf
       + 832/9*theta( - 3 + N)*den( - 1 + N)*S(R(-2), - 1 + N)*cf*ca*nf
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(-2,1), - 1 + N)*cf*ca*nf
       - 2284/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*ca*nf
       - 160/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf*nf^2
       + 2092/27*theta( - 3 + N)*den( - 1 + N)*S(R(1), - 1 + N)*cf^2*nf
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,-2), - 1 + N)*cf*ca*nf
       + 16/9*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*ca*nf
       - 32/9*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf*nf^2
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1), - 1 + N)*cf^2*nf
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,1,1), - 1 + N)*cf^2*nf
       + 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf*ca*nf
       - 64/3*theta( - 3 + N)*den( - 1 + N)*S(R(1,2), - 1 + N)*cf^2*nf
       + 24*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf*ca*nf
       + 368/9*theta( - 3 + N)*den( - 1 + N)*S(R(2), - 1 + N)*cf^2*nf
       + 32*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf*ca*nf
       - 32/3*theta( - 3 + N)*den( - 1 + N)*S(R(3), - 1 + N)*cf^2*nf
       + 48*theta( - 3 + N)*den(1 + N)*cf*ca*nf*z3
       + 2528/3*theta( - 3 + N)*den(1 + N)*cf*ca*nf
       + 1480/27*theta( - 3 + N)*den(1 + N)*cf*nf^2
       - 48*theta( - 3 + N)*den(1 + N)*cf^2*nf*z3
       - 3338/27*theta( - 3 + N)*den(1 + N)*cf^2*nf
       + 96*theta( - 3 + N)*den(1 + N)^2*cf*ca*nf*z3
       - 13216/27*theta( - 3 + N)*den(1 + N)^2*cf*ca*nf
       - 560/27*theta( - 3 + N)*den(1 + N)^2*cf*nf^2
       - 96*theta( - 3 + N)*den(1 + N)^2*cf^2*nf*z3
       + 5486/9*theta( - 3 + N)*den(1 + N)^2*cf^2*nf
       + 1588/9*theta( - 3 + N)*den(1 + N)^3*cf*ca*nf
       - 232/9*theta( - 3 + N)*den(1 + N)^3*cf*nf^2
       - 212/3*theta( - 3 + N)*den(1 + N)^3*cf^2*nf
       + 248/3*theta( - 3 + N)*den(1 + N)^4*cf*ca*nf
       + 16/3*theta( - 3 + N)*den(1 + N)^4*cf*nf^2
       + 40*theta( - 3 + N)*den(1 + N)^4*cf^2*nf
       + 48*theta( - 3 + N)*den(1 + N)^5*cf*ca*nf
       + 64*theta( - 3 + N)*den(1 + N)^5*cf^2*nf
       - 8*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf*ca*nf
       - 96*theta( - 3 + N)*den(1 + N)^4*S(R(1),1 + N)*cf^2*nf
       - 16*theta( - 3 + N)*den(1 + N)^3*S(R(-2),1 + N)*cf*ca*nf
       - 284/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf*nf^2
       - 48*theta( - 3 + N)*den(1 + N)^3*S(R(1),1 + N)*cf^2*nf
       - 16*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(1 + N)^3*S(R(1,1),1 + N)*cf^2*nf
       - 64*theta( - 3 + N)*den(1 + N)^3*S(R(2),1 + N)*cf^2*nf
       - 48*theta( - 3 + N)*den(1 + N)^2*S(R(-3),1 + N)*cf*ca*nf
       - 104*theta( - 3 + N)*den(1 + N)^2*S(R(-2),1 + N)*cf*ca*nf
       + 728/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*ca*nf
       + 64/9*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf*nf^2
       - 20*theta( - 3 + N)*den(1 + N)^2*S(R(1),1 + N)*cf^2*nf
       - 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,-2),1 + N)*cf*ca*nf
       - 136/3*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf*ca*nf
       + 16/3*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf*nf^2
       + 64*theta( - 3 + N)*den(1 + N)^2*S(R(1,1),1 + N)*cf^2*nf
       + 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,1,1),1 + N)*cf^2*nf
       - 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(1 + N)^2*S(R(1,2),1 + N)*cf^2*nf
       - 40*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf*ca*nf
       - 112*theta( - 3 + N)*den(1 + N)^2*S(R(2),1 + N)*cf^2*nf
       - 56*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(1 + N)^2*S(R(3),1 + N)*cf^2*nf
       - 8*theta( - 3 + N)*den(1 + N)*S(R(-3),1 + N)*cf*ca*nf
       + 544/3*theta( - 3 + N)*den(1 + N)*S(R(-2),1 + N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(1 + N)*S(R(-2,1),1 + N)*cf*ca*nf
       + 2912/9*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*ca*nf
       - 104/3*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf*nf^2
       - 1004/3*theta( - 3 + N)*den(1 + N)*S(R(1),1 + N)*cf^2*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1,-2),1 + N)*cf*ca*nf
       - 12*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*ca*nf
       + 8/3*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf*nf^2
       + 28/3*theta( - 3 + N)*den(1 + N)*S(R(1,1),1 + N)*cf^2*nf
       + 16*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf*ca*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1,1,1),1 + N)*cf^2*nf
       - 16*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(1 + N)*S(R(1,2),1 + N)*cf^2*nf
       + 40/3*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf*ca*nf
       + 272/3*theta( - 3 + N)*den(1 + N)*S(R(2),1 + N)*cf^2*nf
       - 36*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf*ca*nf
       + 8*theta( - 3 + N)*den(1 + N)*S(R(3),1 + N)*cf^2*nf
       + 64*theta( - 3 + N)*den(2 + N)*cf*ca*nf*z3
       - 82232/81*theta( - 3 + N)*den(2 + N)*cf*ca*nf
       - 800/27*theta( - 3 + N)*den(2 + N)*cf*nf^2
       - 64*theta( - 3 + N)*den(2 + N)*cf^2*nf*z3
       + 872/3*theta( - 3 + N)*den(2 + N)*cf^2*nf
       - 18392/27*theta( - 3 + N)*den(2 + N)^2*cf*ca*nf
       - 304/9*theta( - 3 + N)*den(2 + N)^2*cf*nf^2
       - 104/27*theta( - 3 + N)*den(2 + N)^2*cf^2*nf
       - 664/3*theta( - 3 + N)*den(2 + N)^3*cf*ca*nf
       + 352/3*theta( - 3 + N)*den(2 + N)^3*cf^2*nf
       + 352/3*theta( - 3 + N)*den(2 + N)^4*cf^2*nf
       - 32/3*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*cf*ca*nf
       - 256/3*theta( - 3 + N)*den(2 + N)^3*S(R(1),2 + N)*cf^2*nf
       - 64*theta( - 3 + N)*den(2 + N)^2*S(R(-2),2 + N)*cf*ca*nf
       - 448/9*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf*nf^2
       - 448/9*theta( - 3 + N)*den(2 + N)^2*S(R(1),2 + N)*cf^2*nf
       - 32/3*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den(2 + N)^2*S(R(1,1),2 + N)*cf^2*nf
       - 32/3*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*cf*ca*nf
       - 64*theta( - 3 + N)*den(2 + N)^2*S(R(2),2 + N)*cf^2*nf
       - 128/3*theta( - 3 + N)*den(2 + N)*S(R(-3),2 + N)*cf*ca*nf
       - 1696/9*theta( - 3 + N)*den(2 + N)*S(R(-2),2 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(-2,1),2 + N)*cf*ca*nf
       - 1208/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*ca*nf
       + 592/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf*nf^2
       + 968/27*theta( - 3 + N)*den(2 + N)*S(R(1),2 + N)*cf^2*nf
       - 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,-2),2 + N)*cf*ca*nf
       - 16/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*ca*nf
       + 32/9*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf*nf^2
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,1),2 + N)*cf^2*nf
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf*ca*nf
       - 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,1,1),2 + N)*cf^2*nf
       - 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,2),2 + N)*cf*ca*nf
       + 64/3*theta( - 3 + N)*den(2 + N)*S(R(1,2),2 + N)*cf^2*nf
       - 72*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf*ca*nf
       - 800/9*theta( - 3 + N)*den(2 + N)*S(R(2),2 + N)*cf^2*nf
       - 32*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(2 + N)*S(R(3),2 + N)*cf^2*nf
       - 48*theta( - 3 + N)*den(N)*cf*ca*nf*z3
       - 484/3*theta( - 3 + N)*den(N)*cf*ca*nf
       - 616/27*theta( - 3 + N)*den(N)*cf*nf^2
       + 48*theta( - 3 + N)*den(N)*cf^2*nf*z3
       - 2530/27*theta( - 3 + N)*den(N)*cf^2*nf
       + 96*theta( - 3 + N)*den(N)^2*cf*ca*nf*z3
       - 14260/27*theta( - 3 + N)*den(N)^2*cf*ca*nf
       + 592/27*theta( - 3 + N)*den(N)^2*cf*nf^2
       - 96*theta( - 3 + N)*den(N)^2*cf^2*nf*z3
       + 1510/9*theta( - 3 + N)*den(N)^2*cf^2*nf
       - 1076/9*theta( - 3 + N)*den(N)^3*cf*ca*nf
       - 232/9*theta( - 3 + N)*den(N)^3*cf*nf^2
       + 512/3*theta( - 3 + N)*den(N)^3*cf^2*nf
       - 232/3*theta( - 3 + N)*den(N)^4*cf*ca*nf
       + 16/3*theta( - 3 + N)*den(N)^4*cf*nf^2
       + 24*theta( - 3 + N)*den(N)^4*cf^2*nf
       - 64*theta( - 3 + N)*den(N)^5*cf*ca*nf
       + 64*theta( - 3 + N)*den(N)^5*cf^2*nf
       - 8*theta( - 3 + N)*den(N)^4*S(R(1),N)*cf*ca*nf
       - 96*theta( - 3 + N)*den(N)^4*S(R(1),N)*cf^2*nf
       - 80*theta( - 3 + N)*den(N)^3*S(R(-2),N)*cf*ca*nf
       - 236/3*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf*ca*nf
       + 32/3*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf*nf^2
       + 32*theta( - 3 + N)*den(N)^3*S(R(1),N)*cf^2*nf
       - 16*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(N)^3*S(R(1,1),N)*cf^2*nf
       - 64*theta( - 3 + N)*den(N)^3*S(R(2),N)*cf^2*nf
       - 48*theta( - 3 + N)*den(N)^2*S(R(-3),N)*cf*ca*nf
       - 24*theta( - 3 + N)*den(N)^2*S(R(-2),N)*cf*ca*nf
       + 1916/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*ca*nf
       - 80/9*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf*nf^2
       - 272*theta( - 3 + N)*den(N)^2*S(R(1),N)*cf^2*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(1,-2),N)*cf*ca*nf
       - 88/3*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf*ca*nf
       + 16/3*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf*nf^2
       + 48*theta( - 3 + N)*den(N)^2*S(R(1,1),N)*cf^2*nf
       + 32*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(1,1,1),N)*cf^2*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(N)^2*S(R(1,2),N)*cf^2*nf
       - 16*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf*ca*nf
       - 32*theta( - 3 + N)*den(N)^2*S(R(2),N)*cf^2*nf
       - 56*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(N)^2*S(R(3),N)*cf^2*nf
       + 8*theta( - 3 + N)*den(N)*S(R(-3),N)*cf*ca*nf
       - 256/3*theta( - 3 + N)*den(N)*S(R(-2),N)*cf*ca*nf
       + 32*theta( - 3 + N)*den(N)*S(R(-2,1),N)*cf*ca*nf
       - 1748/9*theta( - 3 + N)*den(N)*S(R(1),N)*cf*ca*nf
       + 56/3*theta( - 3 + N)*den(N)*S(R(1),N)*cf*nf^2
       + 664/3*theta( - 3 + N)*den(N)*S(R(1),N)*cf^2*nf
       + 16*theta( - 3 + N)*den(N)*S(R(1,-2),N)*cf*ca*nf
       + 12*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*ca*nf
       - 8/3*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf*nf^2
       - 28/3*theta( - 3 + N)*den(N)*S(R(1,1),N)*cf^2*nf
       - 16*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf*ca*nf
       + 16*theta( - 3 + N)*den(N)*S(R(1,1,1),N)*cf^2*nf
       + 16*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf*ca*nf
       - 16*theta( - 3 + N)*den(N)*S(R(1,2),N)*cf^2*nf
       + 104/3*theta( - 3 + N)*den(N)*S(R(2),N)*cf*ca*nf
       - 128/3*theta( - 3 + N)*den(N)*S(R(2),N)*cf^2*nf
       + 36*theta( - 3 + N)*den(N)*S(R(3),N)*cf*ca*nf
       - 8*theta( - 3 + N)*den(N)*S(R(3),N)*cf^2*nf
       - 64/3*delta( - 2 + N)*cf*ca*nf*z3
       + 2534/243*delta( - 2 + N)*cf*ca*nf
       - 628/243*delta( - 2 + N)*cf*nf^2
       + 64/3*delta( - 2 + N)*cf^2*nf*z3
       - 3682/243*delta( - 2 + N)*cf^2*nf
      ;

*
********************************************************************************
*

L   Pqg0 =

       + nf * (
          + 2
          - 4*x
          + 4*x^2
          );

L   Pgg0 =

       + ca * (
          - 8
          + 4*[1-x]^-1
          + 4*x^-1
          + 4*x
          - 4*x^2
          )

       + delta_(1 - x)*nf * (
          - 2/3
          )

       + delta_(1 - x)*ca * (
          + 11/3
          );

L   Pgq0 =

       + cf * (
          - 4
          + 4*x^-1
          + 2*x
          );

L   Pqg1 =

       + z2*nf*cf * (
          - 8
          + 16*x
          - 16*x^2
          )

       + z2*nf*ca * (
          - 16*x
          )

       + nf*cf * (
          + 28
          - 58*x
          + 40*x^2
          )

       + nf*ca * (
          - 8
          + 80/9*x^-1
          + 100*x
          - 872/9*x^2
          )

       + H(R(-1,0),x)*nf*ca * (
          - 8
          - 16*x
          - 16*x^2
          )

       + H(R(0),x)*nf*cf * (
          + 6
          - 8*x
          + 16*x^2
          )

       + H(R(0),x)*nf*ca * (
          + 4
          + 32*x
          + 176/3*x^2
          )

       + H(R(0,0),x)*nf*cf * (
          + 4
          - 8*x
          + 16*x^2
          )

       + H(R(0,0),x)*nf*ca * (
          - 8
          - 16*x
          )

       + H(R(1),x)*nf*cf * (
          - 16*x
          + 16*x^2
          )

       + H(R(1),x)*nf*ca * (
          + 16*x
          - 16*x^2
          )

       + H(R(1,0),x)*nf*cf * (
          + 8
          - 16*x
          + 16*x^2
          )

       + H(R(1,1),x)*nf*cf * (
          + 8
          - 16*x
          + 16*x^2
          )

       + H(R(1,1),x)*nf*ca * (
          - 8
          + 16*x
          - 16*x^2
          )

       + H(R(2),x)*nf*cf * (
          + 8
          - 16*x
          + 16*x^2
          );

L   Pgg1 =

       + z2*ca^2 * (
          + 32
          - 8*[1-x]^-1
          - 8*[1+x]^-1
          + 16*x^2
          )

       + nf*cf * (
          - 32
          + 8/3*x^-1
          + 16*x
          + 40/3*x^2
          )

       + nf*ca * (
          + 116/9
          - 40/9*[1-x]^-1
          - 92/9*x^-1
          - 76/9*x
          + 92/9*x^2
          )

       + ca^2 * (
          - 50/9
          + 268/9*[1-x]^-1
          - 218/9*x
          )

       + delta_(1 - x)*z3*ca^2 * (
          + 12
          )

       + delta_(1 - x)*nf*cf * (
          - 2
          )

       + delta_(1 - x)*nf*ca * (
          - 8/3
          )

       + delta_(1 - x)*ca^2 * (
          + 32/3
          )

       + H(R(-1,0),x)*ca^2 * (
          + 32
          - 16*[1+x]^-1
          + 16*x^-1
          + 16*x
          + 16*x^2
          )

       + H(R(0),x)*nf*cf * (
          - 12
          - 20*x
          )

       + H(R(0),x)*nf*ca * (
          - 8/3
          - 8/3*x
          )

       + H(R(0),x)*ca^2 * (
          - 100/3
          + 44/3*x
          - 176/3*x^2
          )

       + H(R(0,0),x)*nf*cf * (
          - 8
          - 8*x
          )

       + H(R(0,0),x)*ca^2 * (
          + 8*[1-x]^-1
          + 8*[1+x]^-1
          + 32*x
          - 16*x^2
          )

       + H(R(1,0),x)*ca^2 * (
          - 32
          + 16*[1-x]^-1
          + 16*x^-1
          + 16*x
          - 16*x^2
          )

       + H(R(2),x)*ca^2 * (
          - 32
          + 16*[1-x]^-1
          + 16*x^-1
          + 16*x
          - 16*x^2
          );

L   Pgq1 =

       + z2*cf*ca * (
          + 16
          )

       + nf*cf * (
          + 80/9
          - 80/9*x^-1
          - 64/9*x
          )

       + cf*ca * (
          + 76/9
          + 4*x^-1
          + 148/9*x
          + 176/9*x^2
          )

       + cf^2 * (
          - 10
          - 14*x
          )

       + H(R(-1,0),x)*cf*ca * (
          + 16
          + 16*x^-1
          + 8*x
          )

       + H(R(0),x)*cf*ca * (
          - 48
          - 20*x
          - 32/3*x^2
          )

       + H(R(0),x)*cf^2 * (
          + 8
          + 14*x
          )

       + H(R(0,0),x)*cf*ca * (
          + 16
          + 8*x
          )

       + H(R(0,0),x)*cf^2 * (
          - 8
          + 4*x
          )

       + H(R(1),x)*nf*cf * (
          - 16/3
          + 16/3*x^-1
          + 8/3*x
          )

       + H(R(1),x)*cf*ca * (
          + 88/3
          - 88/3*x^-1
          - 68/3*x
          )

       + H(R(1),x)*cf^2 * (
          - 24
          + 24*x^-1
          + 20*x
          )

       + H(R(1,0),x)*cf*ca * (
          - 16
          + 16*x^-1
          + 8*x
          )

       + H(R(1,1),x)*cf*ca * (
          - 16
          + 16*x^-1
          + 8*x
          )

       + H(R(1,1),x)*cf^2 * (
          + 16
          - 16*x^-1
          - 8*x
          )

       + H(R(2),x)*cf*ca * (
          - 16
          + 16*x^-1
          + 8*x
          );

L   Pqqps1 =

       + nf*cf * (
          - 8
          + 80/9*x^-1
          + 24*x
          - 224/9*x^2
          )

       + H(R(0),x)*nf*cf * (
          + 4
          + 20*x
          + 32/3*x^2
          )

       + H(R(0,0),x)*nf*cf * (
          - 8
          - 8*x
          );

L   Pgg2 =

       + z2*nf*cf*ca * (
          + 2696/9
          + 608/9*x^-1
          - 3964/9*x
          - 160/9*x^2
          )

       + z2*nf*cf^2 * (
          - 64
          + 1012/3*x
          - 64/3*x^2
          )

       + z2*nf*ca^2 * (
          - 1700/9
          + 160/9*[1-x]^-1
          + 160/9*[1+x]^-1
          - 128/3*x^-1
          + 1076/9*x
          - 176/3*x^2
          )

       + z2*nf^2*cf * (
          + 16/9
          + 64/9*x
          - 32/9*x^2
          )

       + z2*nf^2*ca * (
          + 32/9
          + 32/9*x
          )

       + z2*ca^3 * (
          + 976/3
          - 1072/9*[1-x]^-1
          - 1072/9*[1+x]^-1
          - 3112/9*x^-1
          - 2632/9*x
          + 1072/9*x^2
          )

       + z2^2*nf*cf*ca * (
          - 664/5
          - 344/5*x
          )

       + z2^2*nf*cf^2 * (
          + 104
          + 72*x
          )

       + z2^2*nf*ca^2 * (
          + 276/5
          + 156/5*x
          )

       + z2^2*ca^3 * (
          - 272
          - 24/5*[1-x]^-1
          + 88*[1+x]^-1
          - 464/5*x^-1
          - 1152/5*x
          - 416/5*x^2
          )

       + z3*nf*cf*ca * (
          - 464/3
          + 32*[1-x]^-1
          + 224/3*x^-1
          - 1208/3*x
          + 32/3*x^2
          )

       + z3*nf*cf^2 * (
          + 64/3*x^-1
          - 32*x
          - 256/3*x^2
          )

       + z3*nf*ca^2 * (
          + 496/3
          - 128/3*[1-x]^-1
          - 176/3*x^-1
          + 112/3*x
          + 128/3*x^2
          )

       + z3*nf^2*cf * (
          - 16/3
          - 16/3*x
          )

       + z3*ca^3 * (
          + 1064/3
          + 176/3*[1-x]^-1
          - 1144/3*x^-1
          + 704/3*x
          + 1408/3*x^2
          )

       + nf*cf*ca * (
          + 19204/27
          - 110/3*[1-x]^-1
          - 30662/81*x^-1
          + 182/27*x
          - 24526/81*x^2
          )

       + nf*cf^2 * (
          - 850/3
          - 44/3*x^-1
          + 66*x
          + 232*x^2
          )

       + nf*ca^2 * (
          - 2174/9
          - 836/27*[1-x]^-1
          + 19264/81*x^-1
          + 7358/27*x
          - 19264/81*x^2
          )

       + nf^2*cf * (
          + 224/9
          - 1232/81*x^-1
          + 64/9*x
          - 1360/81*x^2
          )

       + nf^2*ca * (
          - 110/27
          - 16/27*[1-x]^-1
          + 472/81*x^-1
          + 14/3*x
          - 472/81*x^2
          )

       + ca^3 * (
          - 54088/27
          + 490/3*[1-x]^-1
          + 146182/81*x^-1
          + 49678/27*x
          - 146182/81*x^2
          )

       + delta_(1 - x)*z2*z3*ca^3 * (
          - 16
          )

       + delta_(1 - x)*z2*nf*ca^2 * (
          - 8/3
          )

       + delta_(1 - x)*z2*ca^3 * (
          + 8/3
          )

       + delta_(1 - x)*z2^2*nf*ca^2 * (
          - 4/3
          )

       + delta_(1 - x)*z2^2*ca^3 * (
          + 22/3
          )

       + delta_(1 - x)*z3*nf*ca^2 * (
          - 80/3
          )

       + delta_(1 - x)*z3*ca^3 * (
          + 536/3
          )

       + delta_(1 - x)*nf*cf*ca * (
          - 241/18
          )

       + delta_(1 - x)*nf*cf^2 * (
          + 1
          )

       + delta_(1 - x)*nf*ca^2 * (
          - 233/18
          )

       + delta_(1 - x)*nf^2*cf * (
          + 11/9
          )

       + delta_(1 - x)*nf^2*ca * (
          + 29/18
          )

       + delta_(1 - x)*ca^3 * (
          + 79/2
          - 80*z5
          )

       + H(R(-3,0),x)*nf*cf*ca * (
          - 128
          )

       + H(R(-3,0),x)*nf*cf^2 * (
          + 128
          )

       + H(R(-3,0),x)*nf*ca^2 * (
          + 80
          - 16*x
          )

       + H(R(-3,0),x)*ca^3 * (
          - 64
          - 64*[1-x]^-1
          - 64*[1+x]^-1
          + 64*x
          + 128*x^2
          )

       + H(R(-2),x)*z2*nf*cf*ca * (
          + 128
          - 128*x
          )

       + H(R(-2),x)*z2*nf*cf^2 * (
          - 64
          + 64*x
          )

       + H(R(-2),x)*z2*nf*ca^2 * (
          - 48
          + 48*x
          )

       + H(R(-2),x)*z2*ca^3 * (
          - 768
          + 96*[1-x]^-1
          + 256*[1+x]^-1
          - 160*x^-1
          - 96*x
          - 352*x^2
          )

       + H(R(-2,-1,0),x)*nf*cf*ca * (
          + 256
          - 256*x
          )

       + H(R(-2,-1,0),x)*nf*cf^2 * (
          - 128
          + 128*x
          )

       + H(R(-2,-1,0),x)*nf*ca^2 * (
          - 96
          + 96*x
          )

       + H(R(-2,-1,0),x)*ca^3 * (
          - 512
          + 64*[1-x]^-1
          + 128*[1+x]^-1
          - 64*x^-1
          + 64*x
          - 192*x^2
          )

       + H(R(-2,0),x)*nf*cf*ca * (
          - 208
          + 64/3*x^-1
          - 112*x
          + 64*x^2
          )

       + H(R(-2,0),x)*nf*cf^2 * (
          + 128
          - 128/3*x^2
          )

       + H(R(-2,0),x)*nf*ca^2 * (
          + 440/3
          - 32/3*[1-x]^-1
          - 32/3*x^-1
          - 24*x
          - 64/3*x^2
          )

       + H(R(-2,0),x)*ca^3 * (
          + 640/3
          + 176/3*[1-x]^-1
          - 176/3*x^-1
          + 160*x
          + 352*x^2
          )

       + H(R(-2,0,0),x)*nf*cf*ca * (
          - 64
          + 64*x
          )

       + H(R(-2,0,0),x)*nf*cf^2 * (
          + 64
          - 64*x
          )

       + H(R(-2,0,0),x)*nf*ca^2 * (
          + 48
          - 48*x
          )

       + H(R(-2,0,0),x)*ca^3 * (
          + 512
          - 64*[1-x]^-1
          - 288*[1+x]^-1
          + 224*x^-1
          + 416*x
          + 352*x^2
          )

       + H(R(-2,2),x)*ca^3 * (
          + 512
          - 64*[1-x]^-1
          - 192*[1+x]^-1
          + 128*x^-1
          + 128*x
          + 256*x^2
          )

       + H(R(-1),x)*z2*nf*cf*ca * (
          + 208
          - 64/3*x^-1
          + 208*x
          - 64/3*x^2
          )

       + H(R(-1),x)*z2*nf*cf^2 * (
          - 64
          + 64/3*x^-1
          - 64*x
          + 64/3*x^2
          )

       + H(R(-1),x)*z2*nf*ca^2 * (
          - 88
          + 16/3*x^-1
          - 88*x
          + 16/3*x^2
          )

       + H(R(-1),x)*z2*ca^3 * (
          - 424
          - 88*x^-1
          - 424*x
          - 88*x^2
          )

       + H(R(-1),x)*z3*ca^3 * (
          - 384
          + 192*[1+x]^-1
          - 192*x^-1
          - 192*x
          - 192*x^2
          )

       + H(R(-1,-2,0),x)*ca^3 * (
          - 256
          + 128*[1+x]^-1
          - 128*x^-1
          - 128*x
          - 128*x^2
          )

       + H(R(-1,-1),x)*z2*ca^3 * (
          + 512
          - 256*[1+x]^-1
          + 256*x^-1
          + 256*x
          + 256*x^2
          )

       + H(R(-1,-1,0),x)*nf*cf*ca * (
          + 288
          - 256/3*x^-1
          + 288*x
          - 256/3*x^2
          )

       + H(R(-1,-1,0),x)*nf*cf^2 * (
          - 128
          + 128/3*x^-1
          - 128*x
          + 128/3*x^2
          )

       + H(R(-1,-1,0),x)*nf*ca^2 * (
          - 112
          + 32*x^-1
          - 112*x
          + 32*x^2
          )

       + H(R(-1,-1,0),x)*ca^3 * (
          - 48
          + 176*x^-1
          - 48*x
          + 176*x^2
          )

       + H(R(-1,-1,0,0),x)*ca^3 * (
          - 768
          + 384*[1+x]^-1
          - 384*x^-1
          - 384*x
          - 384*x^2
          )

       + H(R(-1,-1,2),x)*ca^3 * (
          - 512
          + 256*[1+x]^-1
          - 256*x^-1
          - 256*x
          - 256*x^2
          )

       + H(R(-1,0),x)*z2*ca^3 * (
          - 576
          + 288*[1+x]^-1
          - 288*x^-1
          - 288*x
          - 288*x^2
          )

       + H(R(-1,0),x)*nf*cf*ca * (
          - 1232/3
          + 704/9*x^-1
          - 1520/3*x
          - 160/9*x^2
          )

       + H(R(-1,0),x)*nf*cf^2 * (
          + 1072/3
          + 1072/3*x
          )

       + H(R(-1,0),x)*nf*ca^2 * (
          + 1208/9
          + 320/9*[1+x]^-1
          - 304/3*x^-1
          + 1528/9*x
          - 304/3*x^2
          )

       + H(R(-1,0),x)*ca^3 * (
          - 872/9
          - 2144/9*[1+x]^-1
          - 680/3*x^-1
          - 3016/9*x
          - 680/3*x^2
          )

       + H(R(-1,0,0),x)*nf*cf*ca * (
          - 144
          - 32/3*x^-1
          - 144*x
          - 32/3*x^2
          )

       + H(R(-1,0,0),x)*nf*cf^2 * (
          + 64
          - 64/3*x^-1
          + 64*x
          - 64/3*x^2
          )

       + H(R(-1,0,0),x)*nf*ca^2 * (
          + 72
          - 32/3*x^-1
          + 72*x
          - 32/3*x^2
          )

       + H(R(-1,0,0),x)*ca^3 * (
          + 608
          + 1408/3*x^-1
          + 608*x
          + 1408/3*x^2
          )

       + H(R(-1,0,0,0),x)*ca^3 * (
          + 512
          - 256*[1+x]^-1
          + 256*x^-1
          + 256*x
          + 256*x^2
          )

       + H(R(-1,2),x)*nf*cf*ca * (
          - 64
          - 64/3*x^-1
          - 64*x
          - 64/3*x^2
          )

       + H(R(-1,2),x)*nf*ca^2 * (
          + 32
          + 32/3*x^-1
          + 32*x
          + 32/3*x^2
          )

       + H(R(-1,2),x)*ca^3 * (
          + 400
          + 176*x^-1
          + 400*x
          + 176*x^2
          )

       + H(R(-1,2,0),x)*ca^3 * (
          + 128
          - 64*[1+x]^-1
          + 64*x^-1
          + 64*x
          + 64*x^2
          )

       + H(R(-1,3),x)*ca^3 * (
          + 512
          - 256*[1+x]^-1
          + 256*x^-1
          + 256*x
          + 256*x^2
          )

       + H(R(0),x)*z2*nf*cf*ca * (
          + 64/3
          + 64/3*x^-1
          - 488/3*x
          - 32*x^2
          )

       + H(R(0),x)*z2*nf*cf^2 * (
          + 64
          + 80*x
          )

       + H(R(0),x)*z2*nf*ca^2 * (
          - 4/3
          + 16/3*[1+x]^-1
          - 32/3*x^-1
          + 140/3*x
          + 32/3*x^2
          )

       + H(R(0),x)*z2*nf^2*cf * (
          + 32/3
          + 32/3*x
          )

       + H(R(0),x)*z2*ca^3 * (
          + 1480/3
          - 88/3*[1+x]^-1
          - 176/3*x^-1
          + 736/3*x
          + 1936/3*x^2
          )

       + H(R(0),x)*z3*nf*cf*ca * (
          - 16
          + 144*x
          )

       + H(R(0),x)*z3*nf*cf^2 * (
          + 112
          - 16*x
          )

       + H(R(0),x)*z3*nf*ca^2 * (
          - 96*x
          )

       + H(R(0),x)*z3*ca^3 * (
          - 112*[1-x]^-1
          - 80*[1+x]^-1
          - 32*x^-1
          - 160*x
          + 192*x^2
          )

       + H(R(0),x)*nf*cf*ca * (
          + 3220/27
          - 1376/27*x^-1
          + 13804/27*x
          + 112*x^2
          )

       + H(R(0),x)*nf*cf^2 * (
          - 254
          - 1082/3*x
          + 72*x^2
          )

       + H(R(0),x)*nf*ca^2 * (
          + 7472/27
          - 8/3*[1-x]^-1
          + 1136/27*x^-1
          - 3958/27*x
          + 2612/9*x^2
          )

       + H(R(0),x)*nf^2*cf * (
          - 352/27
          + 1088/27*x
          + 176/27*x^2
          )

       + H(R(0),x)*nf^2*ca * (
          + 152/27
          + 86/27*x
          + 104/27*x^2
          )

       + H(R(0),x)*ca^3 * (
          + 30794/27
          + 8/3*[1-x]^-1
          + 6320/27*x^-1
          + 9374/27*x
          + 10604/9*x^2
          )

       + H(R(0,0),x)*z2*nf*cf*ca * (
          + 96
          + 96*x
          )

       + H(R(0,0),x)*z2*nf*cf^2 * (
          + 32
          + 32*x
          )

       + H(R(0,0),x)*z2*nf*ca^2 * (
          - 24
          - 40*x
          )

       + H(R(0,0),x)*z2*ca^3 * (
          + 96
          - 128*[1-x]^-1
          - 128*[1+x]^-1
          - 352*x
          + 256*x^2
          )

       + H(R(0,0),x)*nf*cf*ca * (
          + 4028/9
          + 4832/9*x
          + 1016/9*x^2
          )

       + H(R(0,0),x)*nf*cf^2 * (
          - 328/3
          - 396*x
          + 224/3*x^2
          )

       + H(R(0,0),x)*nf*ca^2 * (
          - 208/9
          - 160/9*[1-x]^-1
          - 160/9*[1+x]^-1
          - 2420/9*x
          + 496/9*x^2
          )

       + H(R(0,0),x)*nf^2*cf * (
          + 184/9
          + 232/9*x
          + 64/9*x^2
          )

       + H(R(0,0),x)*nf^2*ca * (
          - 16/9
          - 16/9*x
          )

       + H(R(0,0),x)*ca^3 * (
          - 11620/9
          + 1072/9*[1-x]^-1
          + 1072/9*[1+x]^-1
          + 1052/3*x
          - 12136/9*x^2
          )

       + H(R(0,0,0),x)*nf*cf*ca * (
          + 232/3
          - 560/3*x
          )

       + H(R(0,0,0),x)*nf*cf^2 * (
          - 24
          - 88*x
          + 32/3*x^2
          )

       + H(R(0,0,0),x)*nf*ca^2 * (
          + 32
          - 152/3*x
          )

       + H(R(0,0,0),x)*nf^2*cf * (
          - 16/3
          - 16/3*x
          )

       + H(R(0,0,0),x)*ca^3 * (
          - 368
          + 464/3*x
          - 704*x^2
          )

       + H(R(0,0,0,0),x)*nf*cf*ca * (
          + 128
          - 192*x
          )

       + H(R(0,0,0,0),x)*nf*cf^2 * (
          - 32
          - 32*x
          )

       + H(R(0,0,0,0),x)*nf*ca^2 * (
          + 16*x
          )

       + H(R(0,0,0,0),x)*ca^3 * (
          - 256
          + 64*[1-x]^-1
          + 64*[1+x]^-1
          + 448*x
          - 128*x^2
          )

       + H(R(1),x)*z2*nf*cf*ca * (
          + 144
          + 128/3*x^-1
          - 144*x
          - 128/3*x^2
          )

       + H(R(1),x)*z2*nf*cf^2 * (
          - 64
          - 64/3*x^-1
          + 64*x
          + 64/3*x^2
          )

       + H(R(1),x)*z2*nf*ca^2 * (
          - 56
          - 16*x^-1
          + 56*x
          + 16*x^2
          )

       + H(R(1),x)*z2*ca^3 * (
          - 24
          - 88*x^-1
          + 24*x
          + 88*x^2
          )

       + H(R(1),x)*z3*ca^3 * (
          + 192
          - 96*[1-x]^-1
          - 96*x^-1
          - 96*x
          + 96*x^2
          )

       + H(R(1),x)*nf*cf*ca * (
          - 2264/9
          - 2084/27*x^-1
          + 3428/9*x
          - 1408/27*x^2
          )

       + H(R(1),x)*nf*cf^2 * (
          + 88/3
          + 124/3*x^-1
          - 428/3*x
          + 72*x^2
          )

       + H(R(1),x)*nf*ca^2 * (
          + 1354/9
          + 820/27*x^-1
          - 1354/9*x
          - 820/27*x^2
          )

       + H(R(1),x)*nf^2*cf * (
          - 104/3
          + 256/27*x^-1
          + 56/3*x
          + 176/27*x^2
          )

       + H(R(1),x)*nf^2*ca * (
          + 22/9
          - 104/27*x^-1
          - 22/9*x
          + 104/27*x^2
          )

       + H(R(1),x)*ca^3 * (
          + 124/9
          - 1652/27*x^-1
          - 124/9*x
          + 1652/27*x^2
          )

       + H(R(1,-2,0),x)*ca^3 * (
          + 128
          - 64*[1-x]^-1
          - 64*x^-1
          - 64*x
          + 64*x^2
          )

       + H(R(1,0),x)*z2*ca^3 * (
          + 192
          - 96*[1-x]^-1
          - 96*x^-1
          - 96*x
          + 96*x^2
          )

       + H(R(1,0),x)*nf*cf*ca * (
          - 256/3
          + 136/9*x^-1
          + 112/3*x
          + 296/9*x^2
          )

       + H(R(1,0),x)*nf*cf^2 * (
          - 184/3
          - 80/3*x^-1
          + 40/3*x
          + 224/3*x^2
          )

       + H(R(1,0),x)*nf*ca^2 * (
          + 784/9
          - 320/9*[1-x]^-1
          - 176/3*x^-1
          - 464/9*x
          + 176/3*x^2
          )

       + H(R(1,0),x)*nf^2*cf * (
          - 16/3
          - 64/9*x^-1
          + 16/3*x
          + 64/9*x^2
          )

       + H(R(1,0),x)*ca^3 * (
          - 2344/9
          + 2144/9*[1-x]^-1
          + 1072/9*x^-1
          + 200/9*x
          - 1072/9*x^2
          )

       + H(R(1,0,0),x)*nf*cf*ca * (
          - 104
          - 128/3*x^-1
          + 104*x
          + 128/3*x^2
          )

       + H(R(1,0,0),x)*nf*cf^2 * (
          + 56
          + 32/3*x^-1
          - 56*x
          - 32/3*x^2
          )

       + H(R(1,0,0),x)*nf*ca^2 * (
          + 116/3
          - 16/3*[1-x]^-1
          - 100/3*x
          )

       + H(R(1,0,0),x)*ca^3 * (
          - 776/3
          + 88/3*[1-x]^-1
          + 880/3*x^-1
          + 688/3*x
          - 880/3*x^2
          )

       + H(R(1,0,0,0),x)*ca^3 * (
          - 256
          + 128*[1-x]^-1
          + 128*x^-1
          + 128*x
          - 128*x^2
          )

       + H(R(1,1),x)*nf*cf*ca * (
          + 12
          + 16/9*x^-1
          - 12*x
          - 16/9*x^2
          )

       + H(R(1,1),x)*nf*cf^2 * (
          - 28/3
          - 64/3*x^-1
          + 28/3*x
          + 64/3*x^2
          )

       + H(R(1,1),x)*nf^2*cf * (
          - 8/3
          - 32/9*x^-1
          + 8/3*x
          + 32/9*x^2
          )

       + H(R(1,1,0),x)*nf*cf*ca * (
          - 16
          - 64/3*x^-1
          + 16*x
          + 64/3*x^2
          )

       + H(R(1,1,0),x)*nf*cf^2 * (
          + 16
          + 64/3*x^-1
          - 16*x
          - 64/3*x^2
          )

       + H(R(1,1,0,0),x)*ca^3 * (
          - 256
          + 128*[1-x]^-1
          + 128*x^-1
          + 128*x
          - 128*x^2
          )

       + H(R(1,1,1),x)*nf*cf*ca * (
          - 16
          - 64/3*x^-1
          + 16*x
          + 64/3*x^2
          )

       + H(R(1,1,1),x)*nf*cf^2 * (
          + 16
          + 64/3*x^-1
          - 16*x
          - 64/3*x^2
          )

       + H(R(1,2,0),x)*ca^3 * (
          - 256
          + 128*[1-x]^-1
          + 128*x^-1
          + 128*x
          - 128*x^2
          )

       + H(R(1,3),x)*ca^3 * (
          - 256
          + 128*[1-x]^-1
          + 128*x^-1
          + 128*x
          - 128*x^2
          )

       + H(R(2),x)*z2*nf*cf*ca * (
          + 128
          + 128*x
          )

       + H(R(2),x)*z2*nf*cf^2 * (
          - 64
          - 64*x
          )

       + H(R(2),x)*z2*nf*ca^2 * (
          - 48
          - 48*x
          )

       + H(R(2),x)*z2*ca^3 * (
          - 64*[1-x]^-1
          + 32*[1+x]^-1
          - 96*x^-1
          - 160*x
          + 32*x^2
          )

       + H(R(2),x)*nf*cf*ca * (
          - 2696/9
          + 32/3*x^-1
          - 596/9*x
          + 160/9*x^2
          )

       + H(R(2),x)*nf*cf^2 * (
          + 64
          + 20*x
          + 64/3*x^2
          )

       + H(R(2),x)*nf*ca^2 * (
          + 1700/9
          - 320/9*[1-x]^-1
          - 176/3*x^-1
          + 452/9*x
          + 176/3*x^2
          )

       + H(R(2),x)*nf^2*cf * (
          - 16/9
          - 64/9*x
          + 32/9*x^2
          )

       + H(R(2),x)*nf^2*ca * (
          - 32/9
          - 32/9*x
          )

       + H(R(2),x)*ca^3 * (
          - 976/3
          + 2144/9*[1-x]^-1
          + 1072/9*x^-1
          - 128/3*x
          - 1072/9*x^2
          )

       + H(R(2,0),x)*nf*cf*ca * (
          - 16/3
          - 64/3*x^-1
          + 8/3*x
          + 32/3*x^2
          )

       + H(R(2,0),x)*nf*cf^2 * (
          - 80
          - 128*x
          - 64/3*x^2
          )

       + H(R(2,0),x)*nf*ca^2 * (
          - 32/3
          - 32/3*x
          )

       + H(R(2,0),x)*nf^2*cf * (
          - 32/3
          - 32/3*x
          )

       + H(R(2,0),x)*ca^3 * (
          - 688/3
          + 352/3*x^-1
          + 464/3*x
          - 352*x^2
          )

       + H(R(2,0,0),x)*nf*cf*ca * (
          - 112
          - 112*x
          )

       + H(R(2,0,0),x)*nf*cf^2 * (
          + 48
          + 48*x
          )

       + H(R(2,0,0),x)*nf*ca^2 * (
          + 24
          + 24*x
          )

       + H(R(2,0,0),x)*ca^3 * (
          - 160
          + 160*[1-x]^-1
          + 160*x^-1
          + 320*x
          - 160*x^2
          )

       + H(R(2,1),x)*nf*cf*ca * (
          + 88/3
          - 32/3*x^-1
          + 136/3*x
          + 32/3*x^2
          )

       + H(R(2,1),x)*nf*cf^2 * (
          - 48
          - 64*x
          - 64/3*x^2
          )

       + H(R(2,1),x)*nf^2*cf * (
          - 16/3
          - 16/3*x
          )

       + H(R(2,1,0),x)*nf*cf*ca * (
          - 32
          - 32*x
          )

       + H(R(2,1,0),x)*nf*cf^2 * (
          + 32
          + 32*x
          )

       + H(R(2,1,0),x)*ca^3 * (
          - 256
          + 128*[1-x]^-1
          + 128*x^-1
          + 128*x
          - 128*x^2
          )

       + H(R(2,1,1),x)*nf*cf*ca * (
          - 32
          - 32*x
          )

       + H(R(2,1,1),x)*nf*cf^2 * (
          + 32
          + 32*x
          )

       + H(R(2,2),x)*ca^3 * (
          - 256
          + 128*[1-x]^-1
          + 128*x^-1
          + 128*x
          - 128*x^2
          )

       + H(R(3),x)*nf*cf*ca * (
          - 64/3
          + 152/3*x
          + 32*x^2
          )

       + H(R(3),x)*nf*cf^2 * (
          - 64
          - 80*x
          )

       + H(R(3),x)*nf*ca^2 * (
          + 4/3
          - 16/3*[1-x]^-1
          - 212/3*x
          - 32/3*x^2
          )

       + H(R(3),x)*nf^2*cf * (
          - 32/3
          - 32/3*x
          )

       + H(R(3),x)*ca^3 * (
          - 1480/3
          + 88/3*[1-x]^-1
          - 256/3*x
          - 1936/3*x^2
          )

       + H(R(3,0),x)*nf*cf*ca * (
          - 64
          - 64*x
          )

       + H(R(3,0),x)*ca^3 * (
          - 128
          + 160*[1-x]^-1
          + 32*[1+x]^-1
          + 128*x^-1
          + 384*x
          - 192*x^2
          )

       + H(R(3,1),x)*nf*cf*ca * (
          - 16
          - 16*x
          )

       + H(R(3,1),x)*nf*cf^2 * (
          + 16
          + 16*x
          )

       + H(R(3,1),x)*ca^3 * (
          - 256
          + 128*[1-x]^-1
          + 128*x^-1
          + 128*x
          - 128*x^2
          )

       + H(R(4),x)*nf*cf*ca * (
          - 96
          - 96*x
          )

       + H(R(4),x)*nf*cf^2 * (
          - 32
          - 32*x
          )

       + H(R(4),x)*nf*ca^2 * (
          + 24
          + 24*x
          )

       + H(R(4),x)*ca^3 * (
          - 96
          + 128*[1-x]^-1
          + 128*[1+x]^-1
          + 416*x
          - 256*x^2
          );

L   Pqg2 =

       + z2*nf*cf*ca * (
          + 688/9
          + 268/9*x
          + 392*x^2
          )

       + z2*nf*cf^2 * (
          - 68
          + 48*x
          - 244*x^2
          )

       + z2*nf*ca^2 * (
          - 628/3
          + 512/9*x^-1
          + 2096/9*x
          - 620/9*x^2
          )

       + z2*nf^2*cf * (
          + 80/9
          - 112/9*x
          + 112/9*x^2
          )

       + z2*nf^2*ca * (
          + 88/9*x
          )

       + z2^2*nf*cf*ca * (
          - 292/5
          - 808/5*x
          - 32*x^2
          )

       + z2^2*nf*cf^2 * (
          + 32
          - 192/5*x
          + 352/5*x^2
          )

       + z2^2*nf*ca^2 * (
          + 538/5
          + 396*x
          + 224/5*x^2
          )

       + z3*nf*cf*ca * (
          - 100/3
          - 64*x^-1
          - 1792/3*x
          + 2240/3*x^2
          )

       + z3*nf*cf^2 * (
          - 24
          + 136*x
          - 304*x^2
          )

       + z3*nf*ca^2 * (
          + 48
          + 96*x^-1
          + 1000*x
          - 2768/3*x^2
          )

       + z3*nf^2*cf * (
          - 8/3
          + 16/3*x
          - 16/3*x^2
          )

       + z3*nf^2*ca * (
          - 16*x
          )

       + nf*cf*ca * (
          + 1429/12
          + 220/3*x^-1
          + 2303/2*x
          - 3997/3*x^2
          )

       + nf*cf^2 * (
          + 963/4
          - 641/2*x
          + 81*x^2
          )

       + nf*ca^2 * (
          + 10454/27
          - 9404/27*x^-1
          - 102608/27*x
          + 103040/27*x^2
          )

       + nf^2*cf * (
          - 13025/54
          + 1232/81*x^-1
          - 89/27*x
          + 18424/81*x^2
          )

       + nf^2*ca * (
          + 226/27
          - 424/81*x^-1
          - 1348/27*x
          + 2800/81*x^2
          )

       + H(R(-3,0),x)*nf*cf*ca * (
          + 48
          - 32*x
          )

       + H(R(-3,0),x)*nf*cf^2 * (
          - 64
          - 128*x^2
          )

       + H(R(-3,0),x)*nf*ca^2 * (
          + 72
          - 16*x
          )

       + H(R(-2),x)*z2*nf*cf*ca * (
          + 112
          + 128*x
          + 288*x^2
          )

       + H(R(-2),x)*z2*nf*cf^2 * (
          + 32
          )

       + H(R(-2),x)*z2*nf*ca^2 * (
          + 56
          - 16*x
          + 64*x^2
          )

       + H(R(-2,-1,0),x)*nf*cf*ca * (
          + 32
          + 128*x
          + 192*x^2
          )

       + H(R(-2,-1,0),x)*nf*cf^2 * (
          + 64
          )

       + H(R(-2,-1,0),x)*nf*ca^2 * (
          + 48
          - 160*x
          )

       + H(R(-2,0),x)*nf*cf*ca * (
          + 208
          - 64*x
          + 80*x^2
          )

       + H(R(-2,0),x)*nf*cf^2 * (
          - 176
          - 32*x
          - 128*x^2
          )

       + H(R(-2,0),x)*nf*ca^2 * (
          - 72
          + 216*x
          - 400*x^2
          )

       + H(R(-2,0,0),x)*nf*cf*ca * (
          - 64
          - 96*x
          - 192*x^2
          )

       + H(R(-2,0,0),x)*nf*cf^2 * (
          - 64
          - 64*x
          - 128*x^2
          )

       + H(R(-2,0,0),x)*nf*ca^2 * (
          + 24
          - 208*x
          - 32*x^2
          )

       + H(R(-2,2),x)*nf*cf*ca * (
          - 96
          - 64*x
          - 192*x^2
          )

       + H(R(-2,2),x)*nf*ca^2 * (
          - 32
          - 64*x
          - 64*x^2
          )

       + H(R(-1),x)*z2*nf*cf*ca * (
          + 40
          + 80*x
          + 40*x^2
          )

       + H(R(-1),x)*z2*nf*cf^2 * (
          + 64
          + 128*x
          + 64*x^2
          )

       + H(R(-1),x)*z2*nf*ca^2 * (
          + 120
          + 32/3*x^-1
          + 88*x
          - 64/3*x^2
          )

       + H(R(-1),x)*z3*nf*cf*ca * (
          + 136
          + 272*x
          + 272*x^2
          )

       + H(R(-1),x)*z3*nf*ca^2 * (
          - 40
          - 80*x
          - 80*x^2
          )

       + H(R(-1,-2,0),x)*nf*cf*ca * (
          + 48
          + 96*x
          + 96*x^2
          )

       + H(R(-1,-2,0),x)*nf*cf^2 * (
          + 32
          + 64*x
          + 64*x^2
          )

       + H(R(-1,-2,0),x)*nf*ca^2 * (
          - 16
          - 32*x
          - 32*x^2
          )

       + H(R(-1,-1),x)*z2*nf*cf*ca * (
          - 144
          - 288*x
          - 288*x^2
          )

       + H(R(-1,-1),x)*z2*nf*ca^2 * (
          + 16
          + 32*x
          + 32*x^2
          )

       + H(R(-1,-1,-1,0),x)*nf*cf*ca * (
          - 96
          - 192*x
          - 192*x^2
          )

       + H(R(-1,-1,-1,0),x)*nf*ca^2 * (
          + 96
          + 192*x
          + 192*x^2
          )

       + H(R(-1,-1,0),x)*nf*cf*ca * (
          - 144
          - 224*x
          - 80*x^2
          )

       + H(R(-1,-1,0),x)*nf*cf^2 * (
          + 128
          + 256*x
          + 128*x^2
          )

       + H(R(-1,-1,0),x)*nf*ca^2 * (
          + 112
          - 64/3*x^-1
          - 80*x
          - 640/3*x^2
          )

       + H(R(-1,-1,0,0),x)*nf*cf*ca * (
          + 96
          + 192*x
          + 192*x^2
          )

       + H(R(-1,-1,0,0),x)*nf*cf^2 * (
          + 64
          + 128*x
          + 128*x^2
          )

       + H(R(-1,-1,0,0),x)*nf*ca^2 * (
          + 32
          + 64*x
          + 64*x^2
          )

       + H(R(-1,-1,2),x)*nf*cf*ca * (
          + 96
          + 192*x
          + 192*x^2
          )

       + H(R(-1,-1,2),x)*nf*ca^2 * (
          + 32
          + 64*x
          + 64*x^2
          )

       + H(R(-1,0),x)*z2*nf*cf*ca * (
          + 144
          + 288*x
          + 288*x^2
          )

       + H(R(-1,0),x)*nf*cf*ca * (
          + 296
          + 152*x
          - 144*x^2
          )

       + H(R(-1,0),x)*nf*cf^2 * (
          - 320
          - 368*x
          - 80*x^2
          )

       + H(R(-1,0),x)*nf*ca^2 * (
          - 704/9
          + 832/9*x^-1
          + 5864/9*x
          + 2216/3*x^2
          )

       + H(R(-1,0),x)*nf^2*ca * (
          + 80/9
          + 112/9*x
          + 112/9*x^2
          )

       + H(R(-1,0,0),x)*nf*cf*ca * (
          + 40
          + 64*x
          + 48*x^2
          )

       + H(R(-1,0,0),x)*nf*cf^2 * (
          - 128
          - 288*x
          - 160*x^2
          )

       + H(R(-1,0,0),x)*nf*ca^2 * (
          - 256/3
          - 128/3*x^-1
          - 1016/3*x
          - 976/3*x^2
          )

       + H(R(-1,0,0),x)*nf^2*ca * (
          + 16/3
          + 32/3*x
          + 32/3*x^2
          )

       + H(R(-1,0,0,0),x)*nf*cf*ca * (
          - 88
          - 176*x
          - 176*x^2
          )

       + H(R(-1,0,0,0),x)*nf*cf^2 * (
          - 16
          - 32*x
          - 32*x^2
          )

       + H(R(-1,0,0,0),x)*nf*ca^2 * (
          - 24
          - 48*x
          - 48*x^2
          )

       + H(R(-1,2),x)*nf*cf*ca * (
          - 112
          - 192*x
          - 80*x^2
          )

       + H(R(-1,2),x)*nf*ca^2 * (
          - 64
          - 64/3*x^-1
          - 128*x
          - 256/3*x^2
          )

       + H(R(-1,2,0),x)*nf*cf*ca * (
          - 32
          - 64*x
          - 64*x^2
          )

       + H(R(-1,2,1),x)*nf*cf*ca * (
          - 32
          - 64*x
          - 64*x^2
          )

       + H(R(-1,2,1),x)*nf*ca^2 * (
          + 32
          + 64*x
          + 64*x^2
          )

       + H(R(-1,3),x)*nf*cf*ca * (
          - 96
          - 192*x
          - 192*x^2
          )

       + H(R(-1,3),x)*nf*ca^2 * (
          - 32
          - 64*x
          - 64*x^2
          )

       + H(R(0),x)*z2*nf*cf*ca * (
          - 40
          - 360*x
          - 1384/3*x^2
          )

       + H(R(0),x)*z2*nf*cf^2 * (
          - 36
          + 128*x
          - 176*x^2
          )

       + H(R(0),x)*z2*nf*ca^2 * (
          - 44
          + 108*x
          - 248/3*x^2
          )

       + H(R(0),x)*z3*nf*cf*ca * (
          - 88
          - 656*x
          + 32*x^2
          )

       + H(R(0),x)*z3*nf*cf^2 * (
          - 56
          + 208*x
          - 224*x^2
          )

       + H(R(0),x)*z3*nf*ca^2 * (
          + 192
          + 576*x
          )

       + H(R(0),x)*nf*cf*ca * (
          + 2696/27
          + 24569/27*x
          + 11332/27*x^2
          )

       + H(R(0),x)*nf*cf^2 * (
          + 170
          + 57*x
          + 348*x^2
          )

       + H(R(0),x)*nf*ca^2 * (
          - 8696/27
          - 896/27*x^-1
          - 30962/27*x
          - 3040*x^2
          )

       + H(R(0),x)*nf^2*cf * (
          - 4856/27
          - 5270/27*x
          - 560/9*x^2
          )

       + H(R(0),x)*nf^2*ca * (
          - 160/27
          + 560/27*x
          - 1832/27*x^2
          )

       + H(R(0,0),x)*z2*nf*cf*ca * (
          + 96
          + 96*x
          + 64*x^2
          )

       + H(R(0,0),x)*z2*nf*cf^2 * (
          - 56
          + 112*x
          - 224*x^2
          )

       + H(R(0,0),x)*z2*nf*ca^2 * (
          - 12
          + 56*x
          - 96*x^2
          )

       + H(R(0,0),x)*nf*cf*ca * (
          - 1267/9
          - 3190/9*x
          - 2204/9*x^2
          )

       + H(R(0,0),x)*nf*cf^2 * (
          + 89
          + 126*x
          + 352*x^2
          )

       + H(R(0,0),x)*nf*ca^2 * (
          + 2396/9
          + 200*x
          + 10244/9*x^2
          )

       + H(R(0,0),x)*nf^2*cf * (
          - 1370/9
          - 884/9*x
          - 352/9*x^2
          )

       + H(R(0,0),x)*nf^2*ca * (
          + 64/9
          + 104/3*x
          + 64/9*x^2
          )

       + H(R(0,0,0),x)*nf*cf*ca * (
          + 44/3
          + 344/3*x
          + 1264/3*x^2
          )

       + H(R(0,0,0),x)*nf*cf^2 * (
          + 12
          + 64*x
          + 240*x^2
          )

       + H(R(0,0,0),x)*nf*ca^2 * (
          - 56/3
          + 56/3*x
          + 96*x^2
          )

       + H(R(0,0,0),x)*nf^2*cf * (
          - 152/3
          + 160/3*x
          + 32*x^2
          )

       + H(R(0,0,0),x)*nf^2*ca * (
          - 16/3
          + 64/3*x
          )

       + H(R(0,0,0,0),x)*nf*cf*ca * (
          - 32
          - 64*x
          )

       + H(R(0,0,0,0),x)*nf*cf^2 * (
          + 16
          + 128*x^2
          )

       + H(R(0,0,0,0),x)*nf*ca^2 * (
          + 64
          - 240*x
          )

       + H(R(0,0,0,0),x)*nf^2*cf * (
          - 32
          + 64*x
          )

       + H(R(1),x)*z2*nf*cf*ca * (
          - 52
          - 40*x
          + 80*x^2
          )

       + H(R(1),x)*z2*nf*cf^2 * (
          + 32
          + 80*x
          - 112*x^2
          )

       + H(R(1),x)*z2*nf*ca^2 * (
          + 212/3
          + 32/3*x^-1
          - 64/3*x
          - 136/3*x^2
          )

       + H(R(1),x)*z2*nf^2*ca * (
          - 8/3
          + 16/3*x
          - 16/3*x^2
          )

       + H(R(1),x)*z3*nf*cf*ca * (
          + 312
          - 624*x
          + 624*x^2
          )

       + H(R(1),x)*z3*nf*cf^2 * (
          - 112
          + 224*x
          - 224*x^2
          )

       + H(R(1),x)*z3*nf*ca^2 * (
          - 248
          + 496*x
          - 496*x^2
          )

       + H(R(1),x)*nf*cf*ca * (
          - 3692/27
          - 1484/27*x^-1
          + 12892/27*x
          - 9620/27*x^2
          )

       + H(R(1),x)*nf*cf^2 * (
          - 32
          - 282*x
          + 348*x^2
          )

       + H(R(1),x)*nf*ca^2 * (
          + 8288/27
          + 644/9*x^-1
          - 7750/27*x
          - 1520/27*x^2
          )

       + H(R(1),x)*nf^2*cf * (
          + 224/27
          - 148/27*x
          + 112/27*x^2
          )

       + H(R(1),x)*nf^2*ca * (
          - 152/27
          + 112/27*x
          - 112/27*x^2
          )

       + H(R(1,-2,0),x)*nf*cf*ca * (
          + 48
          - 96*x
          + 96*x^2
          )

       + H(R(1,-2,0),x)*nf*cf^2 * (
          - 32
          + 64*x
          - 64*x^2
          )

       + H(R(1,-2,0),x)*nf*ca^2 * (
          - 48
          + 96*x
          - 96*x^2
          )

       + H(R(1,0),x)*z2*nf*cf*ca * (
          + 144
          - 288*x
          + 288*x^2
          )

       + H(R(1,0),x)*z2*nf*cf^2 * (
          - 112
          + 224*x
          - 224*x^2
          )

       + H(R(1,0),x)*z2*nf*ca^2 * (
          - 80
          + 160*x
          - 160*x^2
          )

       + H(R(1,0),x)*nf*cf*ca * (
          - 116/9
          + 368/9*x^-1
          + 1888/9*x
          - 212*x^2
          )

       + H(R(1,0),x)*nf*cf^2 * (
          + 160
          - 352*x
          + 272*x^2
          )

       + H(R(1,0),x)*nf*ca^2 * (
          + 92/3
          + 24*x^-1
          + 688/3*x
          - 276*x^2
          )

       + H(R(1,0),x)*nf^2*cf * (
          - 40/9
          + 128/9*x
          - 128/9*x^2
          )

       + H(R(1,0,0),x)*nf*cf*ca * (
          + 472/3
          + 32/3*x^-1
          - 488/3*x
          )

       + H(R(1,0,0),x)*nf*cf^2 * (
          - 32
          - 48*x
          + 80*x^2
          )

       + H(R(1,0,0),x)*nf*ca^2 * (
          - 304/3
          - 32*x^-1
          - 148/3*x
          + 592/3*x^2
          )

       + H(R(1,0,0),x)*nf^2*cf * (
          - 16/3
          + 32/3*x
          - 32/3*x^2
          )

       + H(R(1,0,0),x)*nf^2*ca * (
          - 8/3
          + 16/3*x
          - 16/3*x^2
          )

       + H(R(1,0,0,0),x)*nf*cf*ca * (
          - 24
          + 48*x
          - 48*x^2
          )

       + H(R(1,0,0,0),x)*nf*cf^2 * (
          + 48
          - 96*x
          + 96*x^2
          )

       + H(R(1,0,0,0),x)*nf*ca^2 * (
          + 40
          - 80*x
          + 80*x^2
          )

       + H(R(1,1),x)*z2*nf*cf*ca * (
          + 96
          - 192*x
          + 192*x^2
          )

       + H(R(1,1),x)*z2*nf*cf^2 * (
          - 80
          + 160*x
          - 160*x^2
          )

       + H(R(1,1),x)*z2*nf*ca^2 * (
          - 16
          + 32*x
          - 32*x^2
          )

       + H(R(1,1),x)*nf*cf*ca * (
          - 808/9
          + 368/9*x^-1
          + 4568/9*x
          - 4096/9*x^2
          )

       + H(R(1,1),x)*nf*cf^2 * (
          + 108
          - 296*x
          + 244*x^2
          )

       + H(R(1,1),x)*nf*ca^2 * (
          - 164/9
          - 160/9*x^-1
          - 1904/9*x
          + 188*x^2
          )

       + H(R(1,1),x)*nf^2*cf * (
          - 80/9
          + 112/9*x
          - 112/9*x^2
          )

       + H(R(1,1),x)*nf^2*ca * (
          + 80/9
          - 112/9*x
          + 112/9*x^2
          )

       + H(R(1,1,0),x)*nf*cf*ca * (
          + 4
          - 64/3*x^-1
          - 88*x
          + 280/3*x^2
          )

       + H(R(1,1,0),x)*nf*cf^2 * (
          - 32
          - 48*x
          + 80*x^2
          )

       + H(R(1,1,0),x)*nf*ca^2 * (
          + 92/3
          + 64/3*x^-1
          + 392/3*x
          - 168*x^2
          )

       + H(R(1,1,0),x)*nf^2*ca * (
          - 8/3
          + 16/3*x
          - 16/3*x^2
          )

       + H(R(1,1,0,0),x)*nf*cf*ca * (
          + 16
          - 32*x
          + 32*x^2
          )

       + H(R(1,1,0,0),x)*nf*cf^2 * (
          + 64
          - 128*x
          + 128*x^2
          )

       + H(R(1,1,0,0),x)*nf*ca^2 * (
          - 16
          + 32*x
          - 32*x^2
          )

       + H(R(1,1,1),x)*nf*cf*ca * (
          - 56/3
          - 32/3*x^-1
          + 424/3*x
          - 416/3*x^2
          )

       + H(R(1,1,1),x)*nf*cf^2 * (
          - 4
          - 80*x
          + 96*x^2
          )

       + H(R(1,1,1),x)*nf*ca^2 * (
          + 68/3
          + 32/3*x^-1
          - 184/3*x
          + 128/3*x^2
          )

       + H(R(1,1,1),x)*nf^2*cf * (
          + 8/3
          - 16/3*x
          + 16/3*x^2
          )

       + H(R(1,1,1),x)*nf^2*ca * (
          - 8/3
          + 16/3*x
          - 16/3*x^2
          )

       + H(R(1,1,1,0),x)*nf*cf*ca * (
          + 48
          - 96*x
          + 96*x^2
          )

       + H(R(1,1,1,0),x)*nf*cf^2 * (
          + 16
          - 32*x
          + 32*x^2
          )

       + H(R(1,1,1,0),x)*nf*ca^2 * (
          - 64
          + 128*x
          - 128*x^2
          )

       + H(R(1,1,1,1),x)*nf*cf*ca * (
          - 64
          + 128*x
          - 128*x^2
          )

       + H(R(1,1,1,1),x)*nf*cf^2 * (
          + 32
          - 64*x
          + 64*x^2
          )

       + H(R(1,1,1,1),x)*nf*ca^2 * (
          + 32
          - 64*x
          + 64*x^2
          )

       + H(R(1,1,2),x)*nf*cf*ca * (
          - 48
          + 96*x
          - 96*x^2
          )

       + H(R(1,1,2),x)*nf*cf^2 * (
          + 80
          - 160*x
          + 160*x^2
          )

       + H(R(1,1,2),x)*nf*ca^2 * (
          - 32
          + 64*x
          - 64*x^2
          )

       + H(R(1,2),x)*nf*cf*ca * (
          - 20
          + 152*x
          - 120*x^2
          )

       + H(R(1,2),x)*nf*cf^2 * (
          + 32
          - 208*x
          + 176*x^2
          )

       + H(R(1,2),x)*nf*ca^2 * (
          - 44/3
          + 184/3*x
          - 184/3*x^2
          )

       + H(R(1,2),x)*nf^2*ca * (
          + 8/3
          - 16/3*x
          + 16/3*x^2
          )

       + H(R(1,2,0),x)*nf*cf^2 * (
          + 96
          - 192*x
          + 192*x^2
          )

       + H(R(1,2,0),x)*nf*ca^2 * (
          - 32
          + 64*x
          - 64*x^2
          )

       + H(R(1,2,1),x)*nf*cf*ca * (
          - 96
          + 192*x
          - 192*x^2
          )

       + H(R(1,2,1),x)*nf*cf^2 * (
          + 96
          - 192*x
          + 192*x^2
          )

       + H(R(1,3),x)*nf*cf*ca * (
          - 96
          + 192*x
          - 192*x^2
          )

       + H(R(1,3),x)*nf*cf^2 * (
          + 112
          - 224*x
          + 224*x^2
          )

       + H(R(1,3),x)*nf*ca^2 * (
          + 48
          - 96*x
          + 96*x^2
          )

       + H(R(2),x)*z2*nf*cf*ca * (
          + 16
          - 64*x
          + 128*x^2
          )

       + H(R(2),x)*z2*nf*cf^2 * (
          - 32
          + 128*x
          - 160*x^2
          )

       + H(R(2),x)*z2*nf*ca^2 * (
          + 24
          + 80*x
          )

       + H(R(2),x)*nf*cf*ca * (
          - 688/9
          + 1100/9*x
          - 392*x^2
          )

       + H(R(2),x)*nf*cf^2 * (
          + 68
          - 416*x
          + 244*x^2
          )

       + H(R(2),x)*nf*ca^2 * (
          + 628/3
          + 320/9*x^-1
          + 1256/3*x
          + 620/9*x^2
          )

       + H(R(2),x)*nf^2*cf * (
          - 80/9
          + 112/9*x
          - 112/9*x^2
          )

       + H(R(2),x)*nf^2*ca * (
          + 8/3*x
          )

       + H(R(2,0),x)*nf*cf*ca * (
          + 32
          + 72*x
          + 392*x^2
          )

       + H(R(2,0),x)*nf*cf^2 * (
          + 48
          - 80*x
          + 80*x^2
          )

       + H(R(2,0),x)*nf*ca^2 * (
          + 16
          + 96*x
          - 40/3*x^2
          )

       + H(R(2,0,0),x)*nf*cf*ca * (
          + 88
          + 16*x
          + 96*x^2
          )

       + H(R(2,0,0),x)*nf*cf^2 * (
          + 16
          - 96*x
          + 128*x^2
          )

       + H(R(2,0,0),x)*nf*ca^2 * (
          - 84
          - 200*x
          - 64*x^2
          )

       + H(R(2,1),x)*nf*cf*ca * (
          - 44/3
          + 496/3*x
          + 192*x^2
          )

       + H(R(2,1),x)*nf*cf^2 * (
          + 36
          - 104*x
          + 96*x^2
          )

       + H(R(2,1),x)*nf*ca^2 * (
          + 32/3*x^-1
          - 32*x
          - 848/3*x^2
          )

       + H(R(2,1),x)*nf^2*cf * (
          + 8/3
          - 16/3*x
          + 16/3*x^2
          )

       + H(R(2,1,0),x)*nf*cf*ca * (
          - 192*x
          + 96*x^2
          )

       + H(R(2,1,0),x)*nf*cf^2 * (
          + 32
          - 64*x
          + 32*x^2
          )

       + H(R(2,1,0),x)*nf*ca^2 * (
          + 32
          + 128*x
          )

       + H(R(2,1,1),x)*nf*cf*ca * (
          - 56
          + 16*x
          - 64*x^2
          )

       + H(R(2,1,1),x)*nf*cf^2 * (
          + 40
          - 80*x
          + 64*x^2
          )

       + H(R(2,1,1),x)*nf*ca^2 * (
          + 16
          + 64*x
          )

       + H(R(2,2),x)*nf*cf*ca * (
          - 32*x^2
          )

       + H(R(2,2),x)*nf*cf^2 * (
          + 64
          - 128*x
          + 160*x^2
          )

       + H(R(3),x)*nf*cf*ca * (
          + 40
          + 296*x
          + 1384/3*x^2
          )

       + H(R(3),x)*nf*cf^2 * (
          + 36
          - 160*x
          + 176*x^2
          )

       + H(R(3),x)*nf*ca^2 * (
          + 44
          + 108*x
          + 248/3*x^2
          )

       + H(R(3,0),x)*nf*cf*ca * (
          - 64
          - 192*x
          )

       + H(R(3,0),x)*nf*cf^2 * (
          + 64
          - 128*x
          + 192*x^2
          )

       + H(R(3,1),x)*nf*cf*ca * (
          - 64
          - 96*x
          - 64*x^2
          )

       + H(R(3,1),x)*nf*cf^2 * (
          + 80
          - 160*x
          + 192*x^2
          )

       + H(R(3,1),x)*nf*ca^2 * (
          + 48
          + 128*x
          )

       + H(R(4),x)*nf*cf*ca * (
          - 96
          - 128*x
          - 64*x^2
          )

       + H(R(4),x)*nf*cf^2 * (
          + 56
          - 112*x
          + 224*x^2
          )

       + H(R(4),x)*nf*ca^2 * (
          + 12
          - 72*x
          + 96*x^2
          );

L   Pgq2 =

       + z2*nf*cf*ca * (
          - 256/9
          - 112/3*x^-1
          + 308/9*x
          + 32/3*x^2
          )

       + z2*nf*cf^2 * (
          - 128/9
          + 496/9*x^-1
          - 560/9*x
          )

       + z2*cf*ca^2 * (
          + 568/9
          - 872/3*x^-1
          + 568/9*x
          - 1664/9*x^2
          )

       + z2*cf^2*ca * (
          + 2456/9
          - 160*x^-1
          - 5356/9*x
          + 96*x^2
          )

       + z2*cf^3 * (
          - 164
          + 112*x^-1
          + 520*x
          )

       + z2^2*cf*ca^2 * (
          - 772/5
          - 544/5*x^-1
          - 366/5*x
          )

       + z2^2*cf^2*ca * (
          - 784/5
          - 432/5*x^-1
          - 184/5*x
          )

       + z2^2*cf^3 * (
          + 744/5
          + 512/5*x^-1
          + 12*x
          )

       + z3*nf*cf*ca * (
          + 448/3
          - 80*x^-1
          - 184/3*x
          )

       + z3*nf*cf^2 * (
          - 496/3
          + 96*x^-1
          + 248/3*x
          )

       + z3*cf*ca^2 * (
          + 1424/3
          - 336*x^-1
          + 700/3*x
          + 256/3*x^2
          )

       + z3*cf^2*ca * (
          - 1472/3
          + 216*x^-1
          - 884/3*x
          - 32*x^2
          )

       + z3*cf^3 * (
          + 416
          - 240*x^-1
          + 208*x
          )

       + nf*cf*ca * (
          - 2320/9
          + 1934/9*x^-1
          + 148/3*x
          - 1048/9*x^2
          )

       + nf*cf^2 * (
          + 8707/27
          - 3250/9*x^-1
          + 3301/54*x
          - 224/9*x^2
          )

       + nf^2*cf * (
          + 16/9
          - 16/9*x^-1
          + 32/9*x
          )

       + cf*ca^2 * (
          - 12710/9
          + 138305/81*x^-1
          + 1789/3*x
          - 33680/81*x^2
          )

       + cf^2*ca * (
          - 11629/18
          + 163*x^-1
          + 3589/36*x
          + 56*x^2
          )

       + cf^3 * (
          + 475/2
          - 94*x^-1
          - 363/4*x
          )

       + H(R(-3,0),x)*cf*ca^2 * (
          - 144
          + 8*x
          )

       + H(R(-3,0),x)*cf^2*ca * (
          - 224
          + 16*x
          )

       + H(R(-3,0),x)*cf^3 * (
          + 256
          )

       + H(R(-2),x)*z2*cf*ca^2 * (
          - 432
          - 160*x^-1
          - 56*x
          )

       + H(R(-2),x)*z2*cf^2*ca * (
          + 96
          )

       + H(R(-2),x)*z2*cf^3 * (
          - 64
          )

       + H(R(-2,-1,0),x)*cf*ca^2 * (
          - 352
          - 64*x^-1
          + 16*x
          )

       + H(R(-2,-1,0),x)*cf^2*ca * (
          + 192
          )

       + H(R(-2,-1,0),x)*cf^3 * (
          - 128
          )

       + H(R(-2,0),x)*nf*cf*ca * (
          + 160/3
          - 32/3*x^-1
          - 16*x
          )

       + H(R(-2,0),x)*nf*cf^2 * (
          - 64
          + 64/3*x^-1
          + 32*x
          )

       + H(R(-2,0),x)*cf*ca^2 * (
          + 1160/3
          - 176/3*x^-1
          + 16*x
          + 64*x^2
          )

       + H(R(-2,0),x)*cf^2*ca * (
          - 576
          - 16*x
          )

       + H(R(-2,0),x)*cf^3 * (
          + 480
          + 80*x
          )

       + H(R(-2,0,0),x)*cf*ca^2 * (
          + 176
          + 224*x^-1
          + 152*x
          )

       + H(R(-2,0,0),x)*cf^2*ca * (
          - 96
          )

       + H(R(-2,0,0),x)*cf^3 * (
          + 128
          + 32*x
          )

       + H(R(-2,2),x)*cf*ca^2 * (
          + 256
          + 128*x^-1
          + 64*x
          )

       + H(R(-1),x)*z2*nf*cf*ca * (
          - 128/3
          - 128/3*x^-1
          - 64/3*x
          )

       + H(R(-1),x)*z2*cf*ca^2 * (
          - 520/3
          + 104*x^-1
          - 392/3*x
          - 32/3*x^2
          )

       + H(R(-1),x)*z2*cf^2*ca * (
          + 144
          - 48*x^-1
          + 88*x
          )

       + H(R(-1),x)*z2*cf^3 * (
          - 224
          - 96*x^-1
          - 160*x
          )

       + H(R(-1),x)*z3*cf*ca^2 * (
          - 368
          - 368*x^-1
          - 184*x
          )

       + H(R(-1),x)*z3*cf^2*ca * (
          + 176
          + 176*x^-1
          + 88*x
          )

       + H(R(-1,-2,0),x)*cf*ca^2 * (
          - 160
          - 160*x^-1
          - 80*x
          )

       + H(R(-1,-2,0),x)*cf^2*ca * (
          + 96
          + 96*x^-1
          + 48*x
          )

       + H(R(-1,-2,0),x)*cf^3 * (
          - 64
          - 64*x^-1
          - 32*x
          )

       + H(R(-1,-1),x)*z2*cf*ca^2 * (
          + 416
          + 416*x^-1
          + 208*x
          )

       + H(R(-1,-1),x)*z2*cf^2*ca * (
          - 160
          - 160*x^-1
          - 80*x
          )

       + H(R(-1,-1,-1,0),x)*cf*ca^2 * (
          + 192
          + 192*x^-1
          + 96*x
          )

       + H(R(-1,-1,-1,0),x)*cf^2*ca * (
          - 192
          - 192*x^-1
          - 96*x
          )

       + H(R(-1,-1,0),x)*nf*cf*ca * (
          - 128/3
          - 128/3*x^-1
          - 64/3*x
          )

       + H(R(-1,-1,0),x)*cf*ca^2 * (
          + 176/3
          + 304*x^-1
          - 176/3*x
          + 64/3*x^2
          )

       + H(R(-1,-1,0),x)*cf^2*ca * (
          + 480
          + 96*x^-1
          + 304*x
          )

       + H(R(-1,-1,0),x)*cf^3 * (
          - 448
          - 192*x^-1
          - 320*x
          )

       + H(R(-1,-1,0,0),x)*cf*ca^2 * (
          - 384
          - 384*x^-1
          - 192*x
          )

       + H(R(-1,-1,0,0),x)*cf^2*ca * (
          + 128
          + 128*x^-1
          + 64*x
          )

       + H(R(-1,-1,0,0),x)*cf^3 * (
          - 128
          - 128*x^-1
          - 64*x
          )

       + H(R(-1,-1,2),x)*cf*ca^2 * (
          - 320
          - 320*x^-1
          - 160*x
          )

       + H(R(-1,-1,2),x)*cf^2*ca * (
          + 64
          + 64*x^-1
          + 32*x
          )

       + H(R(-1,0),x)*z2*cf*ca^2 * (
          - 320
          - 320*x^-1
          - 160*x
          )

       + H(R(-1,0),x)*z2*cf^2*ca * (
          - 32
          - 32*x^-1
          - 16*x
          )

       + H(R(-1,0),x)*z2*cf^3 * (
          + 64
          + 64*x^-1
          + 32*x
          )

       + H(R(-1,0),x)*nf*cf*ca * (
          - 80/3
          - 656/9*x^-1
          + 8*x
          - 32/9*x^2
          )

       + H(R(-1,0),x)*nf*cf^2 * (
          + 496/9*x^-1
          - 48*x
          + 64/9*x^2
          )

       + H(R(-1,0),x)*cf*ca^2 * (
          - 920/3
          - 872/3*x^-1
          + 212/3*x
          - 160*x^2
          )

       + H(R(-1,0),x)*cf^2*ca * (
          - 608
          - 32*x^-1
          - 700*x
          )

       + H(R(-1,0),x)*cf^3 * (
          + 688
          + 112*x^-1
          + 608*x
          )

       + H(R(-1,0,0),x)*nf*cf*ca * (
          + 64/3
          + 64/3*x^-1
          + 32/3*x
          )

       + H(R(-1,0,0),x)*cf*ca^2 * (
          + 1160/3
          + 784/3*x^-1
          + 232/3*x
          + 128/3*x^2
          )

       + H(R(-1,0,0),x)*cf^2*ca * (
          - 240
          - 48*x^-1
          - 144*x
          )

       + H(R(-1,0,0),x)*cf^3 * (
          + 384
          + 192*x^-1
          + 224*x
          )

       + H(R(-1,0,0,0),x)*cf*ca^2 * (
          + 176
          + 176*x^-1
          + 88*x
          )

       + H(R(-1,0,0,0),x)*cf^2*ca * (
          + 48
          + 48*x^-1
          + 24*x
          )

       + H(R(-1,0,0,0),x)*cf^3 * (
          + 32
          + 32*x^-1
          + 16*x
          )

       + H(R(-1,2),x)*nf*cf*ca * (
          + 64/3
          + 64/3*x^-1
          + 32/3*x
          )

       + H(R(-1,2),x)*cf*ca^2 * (
          + 608/3
          + 48*x^-1
          + 304/3*x
          + 64/3*x^2
          )

       + H(R(-1,2),x)*cf^2*ca * (
          + 96
          + 96*x^-1
          + 64*x
          )

       + H(R(-1,2,0),x)*cf*ca^2 * (
          + 64
          + 64*x^-1
          + 32*x
          )

       + H(R(-1,2,1),x)*cf*ca^2 * (
          + 64
          + 64*x^-1
          + 32*x
          )

       + H(R(-1,2,1),x)*cf^2*ca * (
          - 64
          - 64*x^-1
          - 32*x
          )

       + H(R(-1,3),x)*cf*ca^2 * (
          + 256
          + 256*x^-1
          + 128*x
          )

       + H(R(0),x)*z2*nf*cf*ca * (
          - 64/3
          - 32/3*x^-1
          - 80/3*x
          )

       + H(R(0),x)*z2*nf*cf^2 * (
          + 32/3
          + 64/3*x^-1
          + 80/3*x
          )

       + H(R(0),x)*z2*cf*ca^2 * (
          + 2140/3
          - 176/3*x^-1
          + 1028/3*x
          + 96*x^2
          )

       + H(R(0),x)*z2*cf^2*ca * (
          - 752/3
          - 440/3*x
          )

       + H(R(0),x)*z2*cf^3 * (
          + 64
          + 36*x
          )

       + H(R(0),x)*z3*cf*ca^2 * (
          - 224
          - 32*x^-1
          - 80*x
          )

       + H(R(0),x)*z3*cf^2*ca * (
          - 64
          + 64*x
          )

       + H(R(0),x)*z3*cf^3 * (
          + 192
          - 48*x
          )

       + H(R(0),x)*nf*cf*ca * (
          + 6656/27
          + 1208/27*x^-1
          + 1688/27*x
          + 2000/27*x^2
          )

       + H(R(0),x)*nf*cf^2 * (
          - 5174/27
          - 1520/27*x^-1
          - 3134/27*x
          - 352/27*x^2
          )

       + H(R(0),x)*cf*ca^2 * (
          + 38224/27
          + 6320/27*x^-1
          + 3202/27*x
          + 9344/27*x^2
          )

       + H(R(0),x)*cf^2*ca * (
          - 9925/27
          + 7061/27*x
          + 2680/27*x^2
          )

       + H(R(0),x)*cf^3 * (
          + 91
          - 245*x
          )

       + H(R(0,0),x)*z2*cf*ca^2 * (
          - 104
          - 92*x
          )

       + H(R(0,0),x)*z2*cf^2*ca * (
          - 16*x
          )

       + H(R(0,0),x)*z2*cf^3 * (
          + 48
          - 24*x
          )

       + H(R(0,0),x)*nf*cf*ca * (
          - 1088/9
          - 668/9*x
          - 32/9*x^2
          )

       + H(R(0,0),x)*nf*cf^2 * (
          + 2236/9
          + 1306/9*x
          + 64/9*x^2
          )

       + H(R(0,0),x)*cf*ca^2 * (
          - 5080/9
          + 392/9*x
          - 616/9*x^2
          )

       + H(R(0,0),x)*cf^2*ca * (
          - 2590/9
          + 4631/9*x
          - 448/9*x^2
          )

       + H(R(0,0),x)*cf^3 * (
          - 70
          - 555*x
          )

       + H(R(0,0,0),x)*nf*cf*ca * (
          + 32/3
          )

       + H(R(0,0,0),x)*nf*cf^2 * (
          + 16/3
          - 32/3*x
          )

       + H(R(0,0,0),x)*cf*ca^2 * (
          - 272/3
          - 224*x
          - 224/3*x^2
          )

       + H(R(0,0,0),x)*cf^2*ca * (
          - 760/3
          + 608/3*x
          - 128/3*x^2
          )

       + H(R(0,0,0),x)*cf^3 * (
          + 16
          - 132*x
          )

       + H(R(0,0,0,0),x)*nf*cf^2 * (
          + 64
          - 32*x
          )

       + H(R(0,0,0,0),x)*cf*ca^2 * (
          - 128
          + 120*x
          )

       + H(R(0,0,0,0),x)*cf^2*ca * (
          + 64
          + 64*x
          )

       + H(R(0,0,0,0),x)*cf^3 * (
          - 32
          - 32*x
          )

       + H(R(1),x)*z2*nf*cf*ca * (
          + 16/3
          - 16/3*x^-1
          - 8/3*x
          )

       + H(R(1),x)*z2*nf*cf^2 * (
          + 32/3
          - 32/3*x^-1
          - 16/3*x
          )

       + H(R(1),x)*z2*cf*ca^2 * (
          - 352/3
          - 16/3*x^-1
          + 500/3*x
          + 32/3*x^2
          )

       + H(R(1),x)*z2*cf^2*ca * (
          + 1144/3
          - 472/3*x^-1
          - 788/3*x
          )

       + H(R(1),x)*z2*cf^3 * (
          - 256
          + 96*x^-1
          + 152*x
          )

       + H(R(1),x)*z3*cf*ca^2 * (
          - 112
          + 112*x^-1
          + 56*x
          )

       + H(R(1),x)*z3*cf^2*ca * (
          - 48
          + 48*x^-1
          + 24*x
          )

       + H(R(1),x)*z3*cf^3 * (
          + 256
          - 256*x^-1
          - 128*x
          )

       + H(R(1),x)*nf*cf*ca * (
          - 1184/27
          + 2008/27*x^-1
          + 1720/27*x
          + 400/27*x^2
          )

       + H(R(1),x)*nf*cf^2 * (
          + 8/27
          - 56/27*x^-1
          - 1360/27*x
          )

       + H(R(1),x)*nf^2*cf * (
          + 80/9
          - 80/9*x^-1
          - 64/9*x
          )

       + H(R(1),x)*cf*ca^2 * (
          + 6260/27
          - 13516/27*x^-1
          - 766/27*x
          - 1744/27*x^2
          )

       + H(R(1),x)*cf^2*ca * (
          - 16148/27
          + 13006/27*x^-1
          + 10936/27*x
          + 1384/27*x^2
          )

       + H(R(1),x)*cf^3 * (
          + 316
          - 94*x^-1
          - 258*x
          )

       + H(R(1,-2,0),x)*cf*ca^2 * (
          + 32
          - 32*x^-1
          - 16*x
          )

       + H(R(1,-2,0),x)*cf^2*ca * (
          - 160
          + 160*x^-1
          + 80*x
          )

       + H(R(1,-2,0),x)*cf^3 * (
          + 192
          - 192*x^-1
          - 96*x
          )

       + H(R(1,0),x)*z2*cf*ca^2 * (
          + 64
          - 64*x^-1
          - 32*x
          )

       + H(R(1,0),x)*z2*cf^2*ca * (
          - 32
          + 32*x^-1
          + 16*x
          )

       + H(R(1,0),x)*z2*cf^3 * (
          + 64
          - 64*x^-1
          - 32*x
          )

       + H(R(1,0),x)*nf*cf*ca * (
          + 320/9
          - 320/9*x^-1
          - 280/9*x
          )

       + H(R(1,0),x)*nf*cf^2 * (
          - 32/3*x^-1
          - 8/3*x
          )

       + H(R(1,0),x)*cf*ca^2 * (
          + 1096/9
          - 284/3*x^-1
          + 400/9*x
          + 1160/9*x^2
          )

       + H(R(1,0),x)*cf^2*ca * (
          - 224
          + 764/3*x^-1
          + 80/3*x
          )

       + H(R(1,0),x)*cf^3 * (
          - 16
          - 56*x^-1
          - 24*x
          )

       + H(R(1,0,0),x)*nf*cf*ca * (
          + 32/3
          - 32/3*x^-1
          - 16/3*x
          )

       + H(R(1,0,0),x)*nf*cf^2 * (
          - 16
          + 16*x^-1
          + 8*x
          )

       + H(R(1,0,0),x)*cf*ca^2 * (
          - 140/3
          + 176*x^-1
          - 332/3*x
          - 64/3*x^2
          )

       + H(R(1,0,0),x)*cf^2*ca * (
          - 384
          + 128*x^-1
          + 232*x
          )

       + H(R(1,0,0),x)*cf^3 * (
          + 192
          - 48*x^-1
          - 108*x
          )

       + H(R(1,0,0,0),x)*cf*ca^2 * (
          - 48
          + 48*x^-1
          + 24*x
          )

       + H(R(1,0,0,0),x)*cf^2*ca * (
          - 80
          + 80*x^-1
          + 40*x
          )

       + H(R(1,1),x)*z2*cf*ca^2 * (
          + 32
          - 32*x^-1
          - 16*x
          )

       + H(R(1,1),x)*z2*cf^2*ca * (
          + 64
          - 64*x^-1
          - 32*x
          )

       + H(R(1,1),x)*z2*cf^3 * (
          - 96
          + 96*x^-1
          + 48*x
          )

       + H(R(1,1),x)*nf*cf*ca * (
          + 880/9
          - 104*x^-1
          - 488/9*x
          + 32/9*x^2
          )

       + H(R(1,1),x)*nf*cf^2 * (
          - 712/9
          + 664/9*x^-1
          + 512/9*x
          )

       + H(R(1,1),x)*nf^2*cf * (
          - 16/3
          + 16/3*x^-1
          + 8/3*x
          )

       + H(R(1,1),x)*cf*ca^2 * (
          - 100/9
          + 260/3*x^-1
          + 914/9*x
          + 592/9*x^2
          )

       + H(R(1,1),x)*cf^2*ca * (
          - 884/9
          + 64/9*x^-1
          - 1448/9*x
          - 416/9*x^2
          )

       + H(R(1,1),x)*cf^3 * (
          + 96
          - 92*x^-1
          + 54*x
          )

       + H(R(1,1,0),x)*nf*cf*ca * (
          - 80/3
          + 80/3*x^-1
          + 40/3*x
          )

       + H(R(1,1,0),x)*nf*cf^2 * (
          - 64/3
          + 64/3*x^-1
          + 32/3*x
          )

       + H(R(1,1,0),x)*cf*ca^2 * (
          + 440/3
          - 440/3*x^-1
          - 412/3*x
          )

       + H(R(1,1,0),x)*cf^2*ca * (
          - 8/3
          + 8/3*x^-1
          + 148/3*x
          )

       + H(R(1,1,0),x)*cf^3 * (
          - 96
          + 96*x^-1
          + 64*x
          )

       + H(R(1,1,0,0),x)*cf*ca^2 * (
          - 64
          + 64*x^-1
          + 32*x
          )

       + H(R(1,1,0,0),x)*cf^2*ca * (
          - 160
          + 160*x^-1
          + 80*x
          )

       + H(R(1,1,0,0),x)*cf^3 * (
          + 96
          - 96*x^-1
          - 48*x
          )

       + H(R(1,1,1),x)*nf*cf*ca * (
          - 80/3
          + 80/3*x^-1
          + 40/3*x
          )

       + H(R(1,1,1),x)*nf*cf^2 * (
          + 80/3
          - 80/3*x^-1
          - 40/3*x
          )

       + H(R(1,1,1),x)*cf*ca^2 * (
          + 632/3
          - 688/3*x^-1
          - 340/3*x
          + 32/3*x^2
          )

       + H(R(1,1,1),x)*cf^2*ca * (
          - 1112/3
          + 1120/3*x^-1
          + 712/3*x
          - 32/3*x^2
          )

       + H(R(1,1,1),x)*cf^3 * (
          + 160
          - 144*x^-1
          - 124*x
          )

       + H(R(1,1,1,0),x)*cf*ca^2 * (
          - 128
          + 128*x^-1
          + 64*x
          )

       + H(R(1,1,1,0),x)*cf^2*ca * (
          + 96
          - 96*x^-1
          - 48*x
          )

       + H(R(1,1,1,0),x)*cf^3 * (
          + 32
          - 32*x^-1
          - 16*x
          )

       + H(R(1,1,1,1),x)*cf*ca^2 * (
          - 64
          + 64*x^-1
          + 32*x
          )

       + H(R(1,1,1,1),x)*cf^2*ca * (
          + 128
          - 128*x^-1
          - 64*x
          )

       + H(R(1,1,1,1),x)*cf^3 * (
          - 64
          + 64*x^-1
          + 32*x
          )

       + H(R(1,1,2),x)*cf*ca^2 * (
          - 128
          + 128*x^-1
          + 64*x
          )

       + H(R(1,1,2),x)*cf^2*ca * (
          + 32
          - 32*x^-1
          - 16*x
          )

       + H(R(1,1,2),x)*cf^3 * (
          + 96
          - 96*x^-1
          - 48*x
          )

       + H(R(1,2),x)*nf*cf*ca * (
          - 80/3
          + 80/3*x^-1
          + 40/3*x
          )

       + H(R(1,2),x)*nf*cf^2 * (
          - 32/3
          + 32/3*x^-1
          + 16/3*x
          )

       + H(R(1,2),x)*cf*ca^2 * (
          + 440/3
          - 440/3*x^-1
          - 412/3*x
          )

       + H(R(1,2),x)*cf^2*ca * (
          - 424/3
          + 328/3*x^-1
          + 332/3*x
          )

       + H(R(1,2),x)*cf^3 * (
          + 32
          + 8*x
          )

       + H(R(1,2,0),x)*cf*ca^2 * (
          - 192
          + 192*x^-1
          + 96*x
          )

       + H(R(1,2,0),x)*cf^3 * (
          + 64
          - 64*x^-1
          - 32*x
          )

       + H(R(1,2,1),x)*cf*ca^2 * (
          - 128
          + 128*x^-1
          + 64*x
          )

       + H(R(1,2,1),x)*cf^2*ca * (
          + 64
          - 64*x^-1
          - 32*x
          )

       + H(R(1,2,1),x)*cf^3 * (
          + 64
          - 64*x^-1
          - 32*x
          )

       + H(R(1,3),x)*cf*ca^2 * (
          - 128
          + 128*x^-1
          + 64*x
          )

       + H(R(2),x)*z2*cf*ca^2 * (
          - 80
          - 64*x^-1
          - 56*x
          )

       + H(R(2),x)*z2*cf^2*ca * (
          + 160
          - 32*x^-1
          - 32*x
          )

       + H(R(2),x)*z2*cf^3 * (
          - 96
          + 16*x
          )

       + H(R(2),x)*nf*cf*ca * (
          + 256/9
          - 320/9*x^-1
          - 236/9*x
          - 32/3*x^2
          )

       + H(R(2),x)*nf*cf^2 * (
          + 128/9
          + 128/9*x
          )

       + H(R(2),x)*cf*ca^2 * (
          - 568/9
          + 68/9*x
          + 1664/9*x^2
          )

       + H(R(2),x)*cf^2*ca * (
          - 2456/9
          + 128*x^-1
          - 944/9*x
          - 96*x^2
          )

       + H(R(2),x)*cf^3 * (
          + 164
          + 88*x
          )

       + H(R(2,0),x)*nf*cf*ca * (
          - 32/3
          + 32/3*x^-1
          + 16/3*x
          )

       + H(R(2,0),x)*nf*cf^2 * (
          - 32/3
          + 16/3*x
          )

       + H(R(2,0),x)*cf*ca^2 * (
          - 784/3
          - 32/3*x^-1
          - 592/3*x
          - 224/3*x^2
          )

       + H(R(2,0),x)*cf^2*ca * (
          - 16/3
          + 96*x^-1
          + 224/3*x
          )

       + H(R(2,0),x)*cf^3 * (
          + 48*x
          )

       + H(R(2,0,0),x)*cf*ca^2 * (
          + 104
          + 96*x^-1
          + 84*x
          )

       + H(R(2,0,0),x)*cf^2*ca * (
          - 224
          + 64*x^-1
          + 64*x
          )

       + H(R(2,0,0),x)*cf^3 * (
          + 80
          - 8*x
          )

       + H(R(2,1),x)*nf*cf*ca * (
          - 48
          + 80/3*x^-1
          + 8*x
          )

       + H(R(2,1),x)*nf*cf^2 * (
          - 16/3
          + 8/3*x
          )

       + H(R(2,1),x)*cf*ca^2 * (
          + 104
          - 544/3*x^-1
          - 132*x
          - 128/3*x^2
          )

       + H(R(2,1),x)*cf^2*ca * (
          - 176/3
          + 144*x^-1
          + 448/3*x
          + 32*x^2
          )

       + H(R(2,1),x)*cf^3 * (
          - 16
          - 52*x
          )

       + H(R(2,1,0),x)*cf*ca^2 * (
          - 96
          + 96*x^-1
          + 48*x
          )

       + H(R(2,1,0),x)*cf^2*ca * (
          - 32
          + 32*x^-1
          + 16*x
          )

       + H(R(2,1,1),x)*cf*ca^2 * (
          - 128
          + 64*x^-1
          + 16*x
          )

       + H(R(2,1,1),x)*cf^2*ca * (
          + 112
          - 64*x^-1
          - 8*x
          )

       + H(R(2,1,1),x)*cf^3 * (
          + 16
          - 8*x
          )

       + H(R(2,2),x)*cf*ca^2 * (
          - 96
          + 96*x^-1
          + 48*x
          )

       + H(R(2,2),x)*cf^2*ca * (
          - 64
          + 32*x^-1
          + 32*x
          )

       + H(R(2,2),x)*cf^3 * (
          + 32
          - 16*x
          )

       + H(R(3),x)*nf*cf*ca * (
          + 64/3
          + 32/3*x
          )

       + H(R(3),x)*nf*cf^2 * (
          - 32/3
          + 16/3*x
          )

       + H(R(3),x)*cf*ca^2 * (
          - 2140/3
          - 980/3*x
          - 96*x^2
          )

       + H(R(3),x)*cf^2*ca * (
          + 752/3
          + 392/3*x
          )

       + H(R(3),x)*cf^3 * (
          - 64
          + 44*x
          )

       + H(R(3,0),x)*cf*ca^2 * (
          + 64
          + 128*x^-1
          + 128*x
          )

       + H(R(3,0),x)*cf^2*ca * (
          - 64
          + 32*x
          )

       + H(R(3,1),x)*cf*ca^2 * (
          - 64
          + 128*x^-1
          + 96*x
          )

       + H(R(3,1),x)*cf^2*ca * (
          - 112
          - 8*x
          )

       + H(R(3,1),x)*cf^3 * (
          + 48
          - 24*x
          )

       + H(R(4),x)*cf*ca^2 * (
          + 104
          + 100*x
          )

       + H(R(4),x)*cf^2*ca * (
          + 32*x
          )

       + H(R(4),x)*cf^3 * (
          - 48
          + 24*x
          );

L   Pqqps2 =

       + z2*nf*cf*ca * (
          - 1916/9
          + 512/9*x^-1
          + 904/9*x
          + 448/9*x^2
          )

       + z2*nf*cf^2 * (
          + 272
          + 20*x
          + 448/9*x^2
          )

       + z2*nf^2*cf * (
          + 80/9
          - 64/9*x
          - 32/3*x^2
          )

       + z2^2*nf*cf*ca * (
          + 428/5
          + 468/5*x
          )

       + z2^2*nf*cf^2 * (
          - 296/5
          - 296/5*x
          )

       + z3*nf*cf*ca * (
          - 16/3
          + 96*x^-1
          + 248/3*x
          - 416/3*x^2
          )

       + z3*nf*cf^2 * (
          - 64*x^-1
          + 160*x
          + 256/3*x^2
          )

       + z3*nf^2*cf * (
          + 16/3
          + 16/3*x
          )

       + nf*cf*ca * (
          + 484/3
          - 27044/81*x^-1
          - 2528/3*x
          + 82232/81*x^2
          )

       + nf*cf^2 * (
          + 2530/27
          + 220/3*x^-1
          + 3338/27*x
          - 872/3*x^2
          )

       + nf^2*cf * (
          + 616/27
          + 64/27*x^-1
          - 1480/27*x
          + 800/27*x^2
          )

       + H(R(-3,0),x)*nf*cf*ca * (
          + 80
          - 16*x
          )

       + H(R(-2),x)*z2*nf*cf*ca * (
          + 16
          - 16*x
          )

       + H(R(-2,-1,0),x)*nf*cf*ca * (
          + 32
          - 32*x
          )

       + H(R(-2,0),x)*nf*cf*ca * (
          - 24
          + 104*x
          - 64*x^2
          )

       + H(R(-2,0,0),x)*nf*cf*ca * (
          + 48
          - 48*x
          )

       + H(R(-1),x)*z2*nf*cf*ca * (
          + 40
          + 32/3*x^-1
          + 40*x
          + 32/3*x^2
          )

       + H(R(-1,-1,0),x)*nf*cf*ca * (
          + 16
          - 64/3*x^-1
          + 16*x
          - 64/3*x^2
          )

       + H(R(-1,0),x)*nf*cf*ca * (
          + 256/3
          + 832/9*x^-1
          + 544/3*x
          + 1696/9*x^2
          )

       + H(R(-1,0,0),x)*nf*cf*ca * (
          + 8
          - 128/3*x^-1
          + 8*x
          - 128/3*x^2
          )

       + H(R(-1,2),x)*nf*cf*ca * (
          - 32
          - 64/3*x^-1
          - 32*x
          - 64/3*x^2
          )

       + H(R(0),x)*z2*nf*cf*ca * (
          - 236/3
          + 28/3*x
          - 32/3*x^2
          )

       + H(R(0),x)*z2*nf*cf^2 * (
          + 32
          - 48*x
          - 256/3*x^2
          )

       + H(R(0),x)*z2*nf^2*cf * (
          + 32/3
          + 32/3*x
          )

       + H(R(0),x)*z3*nf*cf*ca * (
          + 208
          + 144*x
          )

       + H(R(0),x)*z3*nf*cf^2 * (
          - 112
          - 112*x
          )

       + H(R(0),x)*nf*cf*ca * (
          - 14260/27
          - 896/27*x^-1
          - 13216/27*x
          - 18392/27*x^2
          )

       + H(R(0),x)*nf*cf^2 * (
          + 1510/9
          + 5486/9*x
          - 104/27*x^2
          )

       + H(R(0),x)*nf^2*cf * (
          + 592/27
          - 560/27*x
          - 304/9*x^2
          )

       + H(R(0,0),x)*z2*nf*cf*ca * (
          + 8
          - 8*x
          )

       + H(R(0,0),x)*z2*nf*cf^2 * (
          + 96
          + 96*x
          )

       + H(R(0,0),x)*nf*cf*ca * (
          + 1076/9
          - 1588/9*x
          + 664/3*x^2
          )

       + H(R(0,0),x)*nf*cf^2 * (
          - 512/3
          + 212/3*x
          - 352/3*x^2
          )

       + H(R(0,0),x)*nf^2*cf * (
          + 232/9
          + 232/9*x
          )

       + H(R(0,0,0),x)*nf*cf*ca * (
          - 232/3
          + 248/3*x
          )

       + H(R(0,0,0),x)*nf*cf^2 * (
          + 24
          + 40*x
          + 352/3*x^2
          )

       + H(R(0,0,0),x)*nf^2*cf * (
          + 16/3
          + 16/3*x
          )

       + H(R(0,0,0,0),x)*nf*cf*ca * (
          + 64
          - 48*x
          )

       + H(R(0,0,0,0),x)*nf*cf^2 * (
          - 64
          - 64*x
          )

       + H(R(1),x)*z2*nf*cf*ca * (
          + 8
          + 32/3*x^-1
          - 8*x
          - 32/3*x^2
          )

       + H(R(1),x)*nf*cf*ca * (
          + 1748/9
          + 2284/27*x^-1
          - 2912/9*x
          + 1208/27*x^2
          )

       + H(R(1),x)*nf*cf^2 * (
          - 664/3
          - 2092/27*x^-1
          + 1004/3*x
          - 968/27*x^2
          )

       + H(R(1),x)*nf^2*cf * (
          - 56/3
          + 160/27*x^-1
          + 104/3*x
          - 592/27*x^2
          )

       + H(R(1,0),x)*nf*cf*ca * (
          + 104/3
          + 24*x^-1
          + 40/3*x
          - 72*x^2
          )

       + H(R(1,0),x)*nf*cf^2 * (
          - 128/3
          + 368/9*x^-1
          + 272/3*x
          - 800/9*x^2
          )

       + H(R(1,0,0),x)*nf*cf*ca * (
          - 36
          - 32*x^-1
          + 36*x
          + 32*x^2
          )

       + H(R(1,0,0),x)*nf*cf^2 * (
          + 8
          + 32/3*x^-1
          - 8*x
          - 32/3*x^2
          )

       + H(R(1,1),x)*nf*cf*ca * (
          - 12
          - 16/9*x^-1
          + 12*x
          + 16/9*x^2
          )

       + H(R(1,1),x)*nf*cf^2 * (
          + 28/3
          + 64/3*x^-1
          - 28/3*x
          - 64/3*x^2
          )

       + H(R(1,1),x)*nf^2*cf * (
          + 8/3
          + 32/9*x^-1
          - 8/3*x
          - 32/9*x^2
          )

       + H(R(1,1,0),x)*nf*cf*ca * (
          + 16
          + 64/3*x^-1
          - 16*x
          - 64/3*x^2
          )

       + H(R(1,1,0),x)*nf*cf^2 * (
          - 16
          - 64/3*x^-1
          + 16*x
          + 64/3*x^2
          )

       + H(R(1,1,1),x)*nf*cf*ca * (
          + 16
          + 64/3*x^-1
          - 16*x
          - 64/3*x^2
          )

       + H(R(1,1,1),x)*nf*cf^2 * (
          - 16
          - 64/3*x^-1
          + 16*x
          + 64/3*x^2
          )

       + H(R(2),x)*z2*nf*cf*ca * (
          + 16
          + 16*x
          )

       + H(R(2),x)*nf*cf*ca * (
          + 1916/9
          + 320/9*x^-1
          + 728/9*x
          - 448/9*x^2
          )

       + H(R(2),x)*nf*cf^2 * (
          - 272
          - 20*x
          - 448/9*x^2
          )

       + H(R(2),x)*nf^2*cf * (
          - 80/9
          + 64/9*x
          + 32/3*x^2
          )

       + H(R(2,0),x)*nf*cf*ca * (
          + 16
          + 40*x
          + 32/3*x^2
          )

       + H(R(2,0),x)*nf*cf^2 * (
          + 32
          + 112*x
          + 64*x^2
          )

       + H(R(2,0,0),x)*nf*cf*ca * (
          - 56
          - 56*x
          )

       + H(R(2,0,0),x)*nf*cf^2 * (
          + 16
          + 16*x
          )

       + H(R(2,1),x)*nf*cf*ca * (
          - 88/3
          + 32/3*x^-1
          - 136/3*x
          - 32/3*x^2
          )

       + H(R(2,1),x)*nf*cf^2 * (
          + 48
          + 64*x
          + 64/3*x^2
          )

       + H(R(2,1),x)*nf^2*cf * (
          + 16/3
          + 16/3*x
          )

       + H(R(2,1,0),x)*nf*cf*ca * (
          + 32
          + 32*x
          )

       + H(R(2,1,0),x)*nf*cf^2 * (
          - 32
          - 32*x
          )

       + H(R(2,1,1),x)*nf*cf*ca * (
          + 32
          + 32*x
          )

       + H(R(2,1,1),x)*nf*cf^2 * (
          - 32
          - 32*x
          )

       + H(R(3),x)*nf*cf*ca * (
          + 236/3
          + 284/3*x
          + 32/3*x^2
          )

       + H(R(3),x)*nf*cf^2 * (
          - 32
          + 48*x
          + 256/3*x^2
          )

       + H(R(3),x)*nf^2*cf * (
          - 32/3
          - 32/3*x
          )

       + H(R(3,0),x)*nf*cf^2 * (
          - 64
          - 64*x
          )

       + H(R(3,1),x)*nf*cf*ca * (
          + 16
          + 16*x
          )

       + H(R(3,1),x)*nf*cf^2 * (
          - 16
          - 16*x
          )

       + H(R(4),x)*nf*cf*ca * (
          - 8
          - 8*x
          )

       + H(R(4),x)*nf*cf^2 * (
          - 96
          - 96*x
          );

*
********************************************************************************
*
