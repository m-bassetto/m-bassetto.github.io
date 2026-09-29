SUBROUTINE fcn(mm,nn,xx,ff)
!****************************************************************************
!
!Purpose: Compute the value of the system of mm equations at point
! xx(nn), producing the result in ff
!This version: no uncertanty, fixed wages and capital rent
!In order to call this subroutine the following global variables
!should be defined:
! (1) current state (o_b,o_pk)
! (2) current csi structure (o_csi) with the valid
!    - Chebyshev coefficients o_coef_tau, o_coef_G, and o_coef_pk,
!      which represent the functions tau(b,pk), G(b,pk), and pk(b,pk)
!    - Chebyshev domain bmin,bmax,pkmin,pkmax
!      which is used by Chebyshev subroutines
!
! Last changed by Vadym Lepetyuk, September 10, 2004
! Copyright by Marco Bassetto, Florin Bidian, Vadym Lepetyuk
! This code can be freely distributed and modified for research purposes only, 
! provided this copyright notice is included in the modified code. 
! Proper credit should be given in all publications arising from
! modifications of this code; this should include a citation of 
! "Politics and Efficiency of Separating Capital and Ordinary Government Budgets"
! by Marco Bassetto with Thomas J. Sargent
!
!****************************************************************************
!DEC$ REAL:8
USE modelsetup
USE chebyshev
IMPLICIT NONE
INTEGER, INTENT(in):: mm,nn
DOUBLE PRECISION, INTENT(in):: xx(nn)
DOUBLE PRECISION, INTENT(out):: ff(mm)

DOUBLE PRECISION cy,l,gam,k,b,tau,G,pk,&
    cy_G,cy_pk,l_G,l_pk,&
    k_G,k_pk,b_G,b_pk,tau_G,tau_pk
DOUBLE PRECISION m1n,moy,mty,m1t
DOUBLE PRECISION tau_n,tau_bPAST,tau_pkPAST,&
    G_n,G_bPAST,G_pkPAST,&
    pk_n,pk_bPAST,pk_pkPAST,&
    co_n,co_kPAST,co_bPAST,co_pkPAST
DOUBLE PRECISION uyc,uyl,uoc
DOUBLE PRECISION uycc,uycl,uocc
DOUBLE PRECISION vG,vG_n,vpk,vpk_n
DOUBLE PRECISION Ek,Eb,Epk
DOUBLE PRECISION rdeltaqk
DOUBLE PRECISION ssr


!CALL writexx(unit2,'fcn:xx',nn,xx)

cy=xx(i__cy);l=xx(i__l);gam=xx(i__gam)
k=xx(i__k);b=xx(i__b);tau=xx(i__tau)
G=xx(i__G);pk=xx(i__pk);cy_G=xx(i__cy_G);cy_pk=xx(i__cy_pk)
l_G=xx(i__l_G);l_pk=xx(i__l_pk);k_G=xx(i__k_G);k_pk=xx(i__k_pk)
b_G=xx(i__b_G);b_pk=xx(i__b_pk);tau_G=xx(i__tau_G);tau_pk=xx(i__tau_pk)

!Defining some constants
m1n=1./(1.+n) !measure of initial old per young
moy=theta*m1n !measure of survived old per young
mty=1+moy !measure of total population per young

!Getting the "future / next" variables and partials
CALL val1_chebyshev(o_csi,b,pk,o_csi%coef(i__tau,:,:),tau_n,tau_bPAST,tau_pkPAST)
CALL val1_chebyshev(o_csi,b,pk,o_csi%coef(i__G,:,:),G_n,G_bPAST,G_pkPAST)
CALL val1_chebyshev(o_csi,b,pk,o_csi%coef(i__pk,:,:),pk_n,pk_bPAST,pk_pkPAST)

rdeltaqk=r+(1.-deltak)*qk
co_n=(rdeltaqk*k+(p+1.)*b)/theta-tau_n
co_kPAST=rdeltaqk/theta
co_bPAST=(p+1.)/theta-tau_bPAST
co_pkPAST=-tau_pkPAST



ff(1)=cy+(p+1.)*o_b*m1n-tau*moy+qk*k+mty*(qG*G+qpk*gam)-w*l

ff(2)=pk-(1.-deltapk)*m1n*o_pk-gam

uyc=u_y_c(cy,l)
uoc=u_o_c(co_n)
ff(3)=qk*uyc-beta*rdeltaqk*uoc

ff(4)=mrscl(cy,l)+w

ff(5)=tau-(1.+alpha*p)*o_b/(1.+theta+n)-(1.-y)*qG*G-(1.-x)*qpk*gam

ff(6)=p*b-(1.+p)*m1n*o_b-mty*(qG*G+qpk*gam-tau)

ff(7)=cy_G-moy*tau_G+qk*k_G+mty*qG-w*l_G

ff(8)=cy_pk-moy*tau_pk+qk*k_pk+mty*qpk-w*l_pk

uycc=u_y_cc(cy,l)
uycl=u_y_cl(cy,l)
uocc=u_o_cc(co_n)
ff(9)=qk*uycc*cy_G+qk*uycl*l_G&
    -beta*rdeltaqk*uocc*(b_G*co_bPAST+k_G*co_kPAST)

ff(10)=qk*uycc*cy_pk+qk*uycl*l_pk&
    -beta*rdeltaqk*uocc*(b_pk*co_bPAST+k_pk*co_kPAST+co_pkPAST)

ff(11)=mrscl_c(cy,l)*cy_G+mrscl_l(cy,l)*l_G

ff(12)=mrscl_c(cy,l)*cy_pk+mrscl_l(cy,l)*l_pk

ff(13)=tau_G-(1.-y)*qG

ff(14)=tau_pk-(1.-x)*qpk

ff(15)=p*b_G-mty*(qG-tau_G)

ff(16)=p*b_pk-mty*(qpk-tau_pk)

uyl=u_y_l(cy,l)
vG=v_G(G)
vG_n=v_G(G_n)
vpk=v_pk(pk)
vpk_n=v_pk(pk_n)
Eb=uoc*co_bPAST+vG_n*G_bPAST+vpk_n*pk_bPAST
Ek=uoc*co_kPAST
ff(17)=uyc*cy_G+uyl*l_G+vG+tbeta*(Eb*b_G+Ek*k_G)

Epk=uoc*co_pkPAST+vG_n*G_pkPAST+vpk_n*pk_pkPAST
ff(18)=uyc*cy_pk+uyl*l_pk+vpk+tbeta*(Eb*b_pk+Ek*k_pk+Epk)


!CALL writeff(unit2,'fcn:ff',mm,ff)

!ssr=SUM(ff**2)
!WRITE(unit1,*) 'fcn:ssr=', ssr

END SUBROUTINE fcn
