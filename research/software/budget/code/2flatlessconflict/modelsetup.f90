MODULE modelsetup
!****************************************************************************
!
! Purpose:
!  Parameters of the model
!  Parameters of the Chebyshev approximation
!
! Last changed by Vadym Lepetyuk, September 6, 2004
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
IMPLICIT NONE

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Variables and parameters of the model
!______________________________________________________________
DOUBLE PRECISION, PARAMETER :: beta=0.96 ! Pure Discount factor
DOUBLE PRECISION, PARAMETER :: theta=0.4 !Probability of surviving (must be >0)
DOUBLE PRECISION, PARAMETER :: tbeta=beta*theta 
DOUBLE PRECISION, PARAMETER :: n=0.02 !Population growth rate
DOUBLE PRECISION, PARAMETER :: alpha=0.0452 !Debt repayment rate
DOUBLE PRECISION :: x !Share of current period public investment paid by bond issue
DOUBLE PRECISION :: y !Share of current period public consumption paid by bond issue

! Technology and shocks
DOUBLE PRECISION, PARAMETER :: deltak=0.06 !Depreciation rates for k
DOUBLE PRECISION, PARAMETER :: deltapk=0.06 !Depreciation rates for pk
DOUBLE PRECISION, PARAMETER :: qk=1. !Relative price of capital
DOUBLE PRECISION, PARAMETER :: qG=1. !Relative price of public consumption
DOUBLE PRECISION, PARAMETER :: qpk=1. !Relative price of public capital 

! Prices
DOUBLE PRECISION, PARAMETER :: w=10. !Wages
DOUBLE PRECISION, PARAMETER :: r=1./.96+0.06-1.0 !Return on Capital
DOUBLE PRECISION, PARAMETER :: p=1./(r/qk-deltak) !Price of console

! Preferences
!u_y(cy,L)=log(cy)+log(1-L)
!u_o(co)=log(co)
!v(G,pk)=phi*log(G)+eta*log(pk)
DOUBLE PRECISION, PARAMETER :: phi=.1327 !Not recalibrated
DOUBLE PRECISION, PARAMETER :: eta=.0202 !Not recalibrated

! Parameters that control adjustment of the size of the box in which the equilibrium is computed
DOUBLE PRECISION, PARAMETER :: adj=0.05 !adjustment factor
DOUBLE PRECISION, PARAMETER :: adjtol=0.02 ! tolerance of values close to bound

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Global definitions
!______________________________________________________________
TYPE csi_type
	INTEGER:: nodes !number of interpolation nodes for b and pk
	DOUBLE PRECISION:: bmin, bmax !range for b
	DOUBLE PRECISION:: pkmin, pkmax !range for pk
	DOUBLE PRECISION, POINTER, DIMENSION(:):: bnodes,pknodes !interpolation nodes DIMENSION(nodes)
	DOUBLE PRECISION, POINTER, DIMENSION(:,:,:)::var !csi DIMENSION(n_var,nodes,nodes)
	!variables related to Chevyshev interpolation
	INTEGER:: order !order of Chebyshev polynomials (order<=nodes-1)
	DOUBLE PRECISION, POINTER, DIMENSION(:)::roots !roots of Chebyshev polynomial DIMENTION(nodes)
	DOUBLE PRECISION, POINTER, DIMENSION(:,:)::T !T(i,j)=Ti(root j) DIMENSION(0:order,1:nodes)
	DOUBLE PRECISION, POINTER, DIMENSION(:)::T2sum !SUM(T(i,:)**2) DIMENSION(0:order)
	DOUBLE PRECISION, POINTER, DIMENSION(:,:,:)::coef !Chebyshev coefficients DIMENSION(n_var,0:order,0:order)
END TYPE csi_type


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Global variables
!______________________________________________________________
INTEGER, PARAMETER:: n_eqn=18 !Number of equations
INTEGER, PARAMETER:: n_var=18 !Number of unknowns
TYPE (csi_type), POINTER:: o_csi !current csi vector
DOUBLE PRECISION o_b, o_pk !current control state for b and pk
INTEGER, PARAMETER :: & !Order of varibles in csi
    i__cy=1,i__l=2,i__gam=3,i__k=4,&
    i__b=5,i__tau=6,i__G=7,i__pk=8,&
    i__cy_G=9,i__cy_pk=10,i__l_G=11,i__l_pk=12,&
    i__k_G=13,i__k_pk=14,i__b_G=15,i__b_pk=16,&
    i__tau_G=17,i__tau_pk=18
INTEGER :: unit0, unit1, unit2 !Logfile units


CONTAINS

! Utility function for young:
! u_y(c,l) = log(c) + log(1-L)

FUNCTION u_y_c(c,l) !partial derivative w.r.t. consumption
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_c
u_y_c=1./c
END FUNCTION u_y_c

FUNCTION u_y_l(c,l) !partial derivative w.r.t. labor
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_l
u_y_l=-1./(1-l)
END FUNCTION u_y_l

FUNCTION u_y_cc(c,l) !2-nd partial derivative w.r.t. consumption
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_cc
u_y_cc=-c**(-2.)
END FUNCTION u_y_cc

FUNCTION u_y_cl(c,l) !cross partial derivative w.r.t. consumption and labor
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_cl
u_y_cl=0.
END FUNCTION u_y_cl

FUNCTION u_y_ll(c,l) !2-nd partial derivative w.r.t. labor
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_ll
u_y_ll=-(1-l)**(-2.)
END FUNCTION u_y_ll

FUNCTION mrscl(c,l) !MRS c/l for young 
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION mrscl
mrscl=-c/(1.-l)
END FUNCTION mrscl

FUNCTION mrscl_c(c,l) !dMRS/dc for young 
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION mrscl_c
mrscl_c=-1./(1.-l)
END FUNCTION mrscl_c

FUNCTION mrscl_l(c,l) !dMRS/dl for young 
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION mrscl_l
mrscl_l=-c/((1.-l)**2)
END FUNCTION mrscl_l

FUNCTION mrscl_cc(c,l) !d2MRS/dc2 for young 
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION mrscl_cc
mrscl_cc=0.
END FUNCTION mrscl_cc

FUNCTION mrscl_cl(c,l) !d2MRS/dcdl for young 
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION mrscl_cl
mrscl_cl=-1./((1.-l)**2)
END FUNCTION mrscl_cl

FUNCTION mrscl_ll(c,l) !dMRS/dl2 for young 
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION mrscl_ll
mrscl_ll=-2*c/((1.-l)**3)
END FUNCTION mrscl_ll

FUNCTION u_y_ccc(c,l) !3-rd partial derivative w.r.t. consumption
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_ccc
u_y_ccc=2.*c**(-3.)
END FUNCTION u_y_ccc

FUNCTION u_y_ccl(c,l) !3-rd partial derivative w.r.t. cons and l
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_ccl
u_y_ccl=0.0
END FUNCTION u_y_ccl

FUNCTION u_y_cll(c,l) !cross partial derivative w.r.t. cons and l
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_cll
u_y_cll=0.0
END FUNCTION u_y_cll

FUNCTION u_y_lll(c,l) !3-nd partial derivative w.r.t. labor
DOUBLE PRECISION, INTENT(in):: c,l
DOUBLE PRECISION u_y_lll
u_y_lll=-2.*(1-l)**(-3.)
END FUNCTION u_y_lll


! Utility function for old:
! u_o(c) = log(c)

FUNCTION u_o_c(c) !1-st derivative
DOUBLE PRECISION, INTENT(in):: c
DOUBLE PRECISION u_o_c
u_o_c=1./c
END FUNCTION u_o_c

FUNCTION u_o_cc(c) !2-nd derivative
DOUBLE PRECISION, INTENT(in):: c
DOUBLE PRECISION u_o_cc
u_o_cc=-c**(-2.)
END FUNCTION u_o_cc

FUNCTION u_o_ccc(c) !3-rd derivative
DOUBLE PRECISION, INTENT(in):: c
DOUBLE PRECISION u_o_ccc
u_o_ccc=2.*c**(-3.)
END FUNCTION u_o_ccc


! Preferences over public goods
! v(G,pk) = phi*log(G) + eta*log(pk)

FUNCTION v_G(G) !1-st derivative w.r.t. G
DOUBLE PRECISION, INTENT(in):: G
DOUBLE PRECISION v_G
v_G=phi/G
END FUNCTION v_G

FUNCTION v_GG(G) !2-st derivative w.r.t. G
DOUBLE PRECISION, INTENT(in):: G
DOUBLE PRECISION v_GG
v_GG=-phi*G**(-2.)
END FUNCTION v_GG

FUNCTION v_pk(pk) !1-st derivative w.r.t. pk
DOUBLE PRECISION, INTENT(in):: pk
DOUBLE PRECISION v_pk
v_pk=eta/pk
END FUNCTION v_pk

FUNCTION v_pkpk(pk) !2-st derivative w.r.t. pk
DOUBLE PRECISION, INTENT(in):: pk
DOUBLE PRECISION v_pkpk
v_pkpk=-eta*pk**(-2.)
END FUNCTION v_pkpk


FUNCTION utility(xx,xxn)
DOUBLE PRECISION, DIMENSION(n_var), INTENT(in):: xx,xxn
DOUBLE PRECISION utility
DOUBLE PRECISION cy,l,G,pk,k,b
DOUBLE PRECISION co_n,G_n,pk_n,tau_n

cy=xx(i__cy);l=xx(i__l);G=xx(i__G);pk=xx(i__pk)
G_n=xxn(i__G);pk_n=xxn(i__pk)

k=xx(i__k);b=xx(i__b);tau_n=xxn(i__tau)
co_n=((r+(1.-deltak))*qk*k+(p+1.)*b)/theta-tau_n

utility = log(cy) + log(1-l) + tbeta*log(co_n) +&
	phi*log(G) + eta*log(pk) +&
	tbeta*phi*log(G_n) + tbeta*eta*log(pk_n)
    
END FUNCTION utility


END MODULE modelsetup
