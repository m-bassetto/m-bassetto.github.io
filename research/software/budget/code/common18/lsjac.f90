SUBROUTINE lsjac(mm,nn,xx,jac,ldfjac)
!****************************************************************************
!
!Purpose: Compute the Jacobian of the system of mm equations at point
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
INTEGER, INTENT(in):: mm,nn,ldfjac
DOUBLE PRECISION, INTENT(in):: xx(nn)
DOUBLE PRECISION, INTENT(out):: jac(ldfjac,nn)

DOUBLE PRECISION cy,co,l,i,gam,k,b,tau,G,pk,&
    cy_G,cy_pk,co_G,co_pk,l_G,l_pk,i_G,i_pk,&
    k_G,k_pk,b_G,b_pk,tau_G,tau_pk
DOUBLE PRECISION m1n,moy,mty,m1t
DOUBLE PRECISION tau_n,tau_bPAST,tau_pkPAST,&
    d2tau_d2bPAST,d2tau_dbPASTdpkPAST,d2tau_d2pkPAST,&
    co_n,co_kPAST,co_bPAST,co_pkPAST,&
    d2co_d2bPAST,d2co_dbPASTdpkPAST,d2co_d2pkPAST,&
    G_n,G_bPAST,G_pkPAST,&
    d2G_d2bPAST,d2G_dbPASTdpkPAST,d2G_d2pkPAST,&
    pk_n,pk_bPAST,pk_pkPAST,&
    d2pk_d2bPAST,d2pk_dbPASTdpkPAST,d2pk_d2pkPAST,&
    vG_n,vpk_n,vGG_n,vpkpk_n
DOUBLE PRECISION uyc,uyl,uycc,uycl,uyll,uyccc,uyccl,uycll
DOUBLE PRECISION uoc,uocc,uoccc
DOUBLE PRECISION E1_k,E2_k,E3_k,E1_b,E2_b,E3_b,E1_pk,E2_pk,E3_pk
DOUBLE PRECISION mrs_cc, mrs_cl, mrs_ll
DOUBLE PRECISION rdeltaqk
DOUBLE PRECISION absjac


!CALL writexx(unit2,'lsjac:xx',nn,xx)

cy=xx(i__cy);l=xx(i__l);gam=xx(i__gam)
k=xx(i__k);b=xx(i__b);tau=xx(i__tau)
G=xx(i__G);pk=xx(i__pk);cy_G=xx(i__cy_G);cy_pk=xx(i__cy_pk)
l_G=xx(i__l_G);l_pk=xx(i__l_pk);k_G=xx(i__k_G);k_pk=xx(i__k_pk)
b_G=xx(i__b_G);b_pk=xx(i__b_pk);tau_G=xx(i__tau_G);tau_pk=xx(i__tau_pk)

!Defining some constants
m1n=1./(1.+n) !measure of initial old per young
moy=theta*m1n !measure of survived old per young
mty=1+moy !measure of total population per young

!Getting the first (pd) and second partial defivatives (pd2)
!for the variables involved co,G,pk
CALL val2_chebyshev(o_csi,b,pk,o_csi%coef(i__tau,:,:), &
    & tau_n,tau_bPAST,tau_pkPAST, &
    & d2tau_d2bPAST,d2tau_dbPASTdpkPAST,d2tau_d2pkPAST)
CALL val2_chebyshev(o_csi,b,pk,o_csi%coef(i__G,:,:), &
    & G_n,G_bPAST,G_pkPAST, &
    & d2G_d2bPAST,d2G_dbPASTdpkPAST,d2G_d2pkPAST)
CALL val2_chebyshev(o_csi,b,pk,o_csi%coef(i__pk,:,:), &
    & pk_n,pk_bPAST,pk_pkPAST, &
    & d2pk_d2bPAST,d2pk_dbPASTdpkPAST,d2pk_d2pkPAST)

rdeltaqk=r+(1.-deltak)*qk
co_n=(rdeltaqk*k+(p+1.)*b)/theta-tau_n
co_kPAST=rdeltaqk/theta
co_bPAST=(p+1.)/theta-tau_bPAST
co_pkPAST=-tau_pkPAST
d2co_d2bPAST=-d2tau_d2bPAST
d2co_dbPASTdpkPAST=-d2tau_dbPASTdpkPAST
d2co_d2pkPAST=-d2tau_d2pkPAST



!Initialize jac
jac=0.  


!Equation 1
jac(1,i__cy)=1.
jac(1,i__l)=-w
jac(1,i__k)=qk
jac(1,i__gam)=mty*qpk
jac(1,i__tau)=-moy
jac(1,i__G)=mty*qG


!Equation 2
jac(2,i__gam)=-1.
jac(2,i__pk)=1.


!Equation 3
uycc=u_y_cc(cy,l)
jac(3,i__cy)=qk*uycc
uycl=u_y_cl(cy,l)
jac(3,i__l)=qk*uycl

uoc=u_o_c(co_n)
uocc=u_o_cc(co_n)
jac(3,i__k)=-beta*rdeltaqk*uocc*co_kPAST 
jac(3,i__b)=-beta*rdeltaqk*uocc*co_bPAST
jac(3,i__pk)=-beta*rdeltaqk*uocc*co_pkPAST


!Equation 4
jac(4,i__cy)=mrscl_c(cy,l)
jac(4,i__l)=mrscl_l(cy,l)


!Equation 5
jac(5,i__gam)=-(1.-x)*qpk
jac(5,i__tau)=1.
jac(5,i__G)=-(1.-y)*qG


!Equation 6
jac(6,i__gam)=-mty*qpk
jac(6,i__b)=p
jac(6,i__tau)=mty
jac(6,i__G)=-mty*qG


!Equation 7
jac(7,i__cy_G)=1.
jac(7,i__l_G)=-w
jac(7,i__k_G)=qk
jac(7,i__tau_G)=-moy


!Equation 8
jac(8,i__cy_pk)=1.
jac(8,i__l_pk)=-w
jac(8,i__k_pk)=qk
jac(8,i__tau_pk)=-moy


!Equation 9
uyccc=u_y_ccc(cy,l)
uyccl=u_y_ccl(cy,l)
uycll=u_y_cll(cy,l)
uoccc=u_o_ccc(co_n)
jac(9,i__cy)=qk*(uyccc*cy_G+uyccl*l_G)
jac(9,i__l)=qk*(uyccl*cy_G+uycll*l_G)

E1_k=uoccc*co_kPAST*co_bPAST !+uocc*d2co_dkPASTdbPAST
E2_k=uoccc*co_kPAST*co_kPAST !+uocc*d2co_d2kPAST
jac(9,i__k)=-beta*rdeltaqk*(b_G*E1_k+k_G*E2_k)

E1_b=uoccc*co_bPAST*co_bPAST+uocc*d2co_d2bPAST
E2_b=uoccc*co_bPAST*co_kPAST !+uocc*d2co_dkPASTdbPAST
jac(9,i__b)=-beta*rdeltaqk*(b_G*E1_b+k_G*E2_b)

E1_pk=uoccc*co_pkPAST*co_bPAST+uocc*d2co_dbPASTdpkPAST
E2_pk=uoccc*co_pkPAST*co_kPAST !+uocc*d2co_dkPASTdpkPAST
jac(9,i__pk)=-beta*rdeltaqk*(b_G*E1_pk+k_G*E2_pk)

jac(9,i__cy_G)=qk*uycc
jac(9,i__l_G)=qk*uycl
jac(9,i__k_G)=jac(3,i__k)
jac(9,i__b_G)=jac(3,i__b)


!Equation 10
jac(10,i__cy)=qk*(uyccc*cy_pk+uyccl*l_pk)
jac(10,i__l)=qk*(uyccl*cy_pk+uycll*l_pk)

E3_k=uoccc*co_kPAST*co_pkPAST !+uocc*d2co_dkPASTdpkPAST
jac(10,i__k)=-beta*rdeltaqk*(b_pk*E1_k+k_pk*E2_k+E3_k)

E3_b=uoccc*co_bPAST*co_pkPAST+uocc*d2co_dbPASTdpkPAST
jac(10,i__b)=-beta*rdeltaqk*(b_pk*E1_b+k_pk*E2_b+E3_b)

E3_pk=uoccc*co_pkPAST*co_pkPAST+uocc*d2co_d2pkPAST
jac(10,i__pk)=-beta*rdeltaqk*(b_pk*E1_pk+k_pk*E2_pk+E3_pk)

jac(10,i__cy_pk)=qk*uycc
jac(10,i__l_pk)=qk*uycl
jac(10,i__k_pk)=jac(3,i__k)
jac(10,i__b_pk)=jac(3,i__b)


!Equation 11
mrs_cc=mrscl_cc(cy,l)
mrs_cl=mrscl_cl(cy,l)
mrs_ll=mrscl_ll(cy,l)
jac(11,i__cy)=mrs_cc*cy_G+mrs_cl*l_G
jac(11,i__l)=mrs_cl*cy_G+mrs_ll*l_G
jac(11,i__cy_G)=jac(4,i__cy)
jac(11,i__l_G)=jac(4,i__l)


!Equation 12
jac(12,i__cy)=mrs_cc*cy_pk+mrs_cl*l_pk
jac(12,i__l)=mrs_cl*cy_pk+mrs_ll*l_pk
jac(12,i__cy_pk)=jac(4,i__cy)
jac(12,i__l_pk)=jac(4,i__l)


!Equation 13
jac(13,i__tau_G)=1.


!Equation 14
jac(14,i__tau_pk)=1.


!Equation 15
jac(15,i__b_G)=p
jac(15,i__tau_G)=mty


!Equation 16
jac(16,i__b_pk)=p
jac(16,i__tau_pk)=mty


!Equation 17
uyll=u_y_ll(cy,l)
jac(17,i__cy)=uycc*cy_G+uycl*l_G
jac(17,i__l)=uycl*cy_G+uyll*l_G

vG_n=v_G(G_n)
vpk_n=v_pk(pk_n)

vGG_n=v_GG(G_n)
vpkpk_n=v_pkpk(pk_n)

!uoc*d2co_dkPASTdbPAST+vGG_n*G_kPAST*G_bPAST+vpk_n*d2pk_dkPASTdbPAST+vG_n*d2G_dkPASTdbPAST+vpkpk_n*pk_kPAST*pk_bPAST
E1_k=uocc*co_bPAST*co_kPAST
E2_k=uocc*co_kPAST*co_kPAST
jac(17,i__k)=tbeta*(b_G*E1_k+k_G*E2_k)

E1_b=uoc*d2co_d2bPAST+uocc*co_bPAST*co_bPAST+vGG_n*G_bPAST*G_bPAST+&
    &vG_n*d2G_d2bPAST+vpkpk_n*pk_bPAST*pk_bPAST+vpk_n*d2pk_d2bPAST
E2_b=uocc*co_kPAST*co_bPAST
jac(17,i__b)=tbeta*(b_G*E1_b+k_G*E2_b)

E1_pk=uoc*d2co_dbPASTdpkPAST+uocc*co_bPAST*co_pkPAST+vGG_n*G_pkPAST*G_bPAST+&
    &vG_n*d2G_dbPASTdpkPAST+vpkpk_n*pk_pkPAST*pk_bPAST+vpk_n*d2pk_dbPASTdpkPAST
E2_pk=uocc*co_kPAST*co_pkPAST
jac(17,i__pk)=tbeta*(b_G*E1_pk+k_G*E2_pk)

jac(17,i__G)=v_GG(G)

uyc=u_y_c(cy,l)
jac(17,i__cy_G)=uyc
uyl=u_y_l(cy,l)
jac(17,i__l_G)=uyl

jac(17,i__k_G)= tbeta*co_kPAST*uoc
jac(17,i__b_G)= tbeta*(co_bPAST*uoc+G_bPAST*vG_n+pk_bPAST*vpk_n)


!Equation 18
jac(18,i__cy)=uycc*cy_pk+uycl*l_pk
jac(18,i__l)=uycl*cy_pk+uyll*l_pk
jac(18,i__cy_pk)=uyc
jac(18,i__l_pk)=uyl
jac(18,i__k_pk)=jac(17,i__k_G)
jac(18,i__b_pk)=jac(17,i__b_G)

E3_k=uocc*co_pkPAST*co_kPAST
jac(18,i__k)=tbeta*(E1_k*b_pk+E2_k*k_pk+E3_k)

E3_b=uoc*d2co_dbPASTdpkPAST+uocc*co_pkPAST*co_bPAST+vGG_n*G_bPAST*G_pkPAST &
    +vG_n*d2G_dbPASTdpkPAST+vpkpk_n*pk_bPAST*pk_pkPAST+vpk_n*d2pk_dbPASTdpkPAST
jac(18,i__b)=tbeta*(E1_b*b_pk+E2_b*k_pk+E3_b)

E3_pk=uoc*d2co_d2pkPAST+uocc*co_pkPAST*co_pkPAST+vGG_n*G_pkPAST*G_pkPAST &
    +vG_n*d2G_d2pkPAST+vpkpk_n*pk_pkPAST*pk_pkPAST+vpk_n*d2pk_d2pkPAST
jac(18,i__pk)=v_pkpk(pk)+tbeta*(E1_pk*b_pk+E2_pk*k_pk+E3_pk)



!absjac=sum(abs(jac))
!WRITE(unit2,*) '-------- lsjac:absjac --------'
!WRITE(unit2,*) 'lsjac:absjac=', absjac

!write(unit1,*) 'absjac',absjac

END SUBROUTINE lsjac
