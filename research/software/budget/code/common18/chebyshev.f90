MODULE chebyshev
!****************************************************************************
!
! Purpose:
!  Subroutines of Chebyshev approximation
! Details:
!  Judd, "Numerical Methods in Economics", 1999, page 238
!
! Last changed by Vadym Lepetyuk, September 20, 2004
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
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER:: PIconst=3.141592653589793238462643383279502884197


CONTAINS
SUBROUTINE coef_chebyshev(csi,val,coef)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Computes Chebyshev coeficients coef(0:order,0:order)
!given the values of function val(1:nodes,1:nodes)
!at the interpolation nodes
!___________________________________________________________ 
IMPLICIT NONE
TYPE (csi_type), intent(in)::csi
DOUBLE PRECISION, dimension(csi%nodes,csi%nodes), intent(in)::val
DOUBLE PRECISION, dimension(0:csi%order,0:csi%order), intent(out)::coef
INTEGER i1,i2,j1,j2

coef=0.0

DO i1=0,csi%order
DO i2=0,csi%order
    DO j1=1,csi%nodes
    DO j2=1,csi%nodes
        coef(i1,i2)=coef(i1,i2)+csi%T(i1,j1)*csi%T(i2,j2)*val(j1,j2)
    END DO
    END DO
    coef(i1,i2)=coef(i1,i2)/csi%T2sum(i1)/csi%T2sum(i2)
END DO
END DO

END SUBROUTINE coef_chebyshev




SUBROUTINE val0_chebyshev(csi,b,pk,coef,val)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Computes Chebyshev approximation of function value
!at the point (b,pk) given the Chebyshev coefficients
!coef(0:order,0:order)
!___________________________________________________________ 
IMPLICIT NONE
TYPE (csi_type), intent(in)::csi
DOUBLE PRECISION, intent(in)::b,pk
DOUBLE PRECISION, dimension(0:csi%order,0:csi%order), intent(in)::coef
DOUBLE PRECISION, intent(out):: val
INTEGER i,i1,i2
DOUBLE PRECISION a1,a2,d1,d2
DOUBLE PRECISION, dimension(0:csi%order)::T1,T2

IF ((b.GT.(csi%bmax)) .OR. (b.LT.(csi%bmin))) &
	STOP 'ERROR: interpolation bounds for b exceeded'

IF ((pk.GT.(csi%pkmax)) .OR. (pk.LT.(csi%pkmin))) &
	STOP 'ERROR: interpolation bounds for pk exceeded'

!The arguments of the polynomial and their derivative, and the 
!derivative coming from arccos
d1=2./(csi%bmax-csi%bmin)
d2=2./(csi%pkmax-csi%pkmin)
a1=(b-csi%bmin)*d1-1.
a2=(pk-csi%pkmin)*d2-1.

!Store once and use multiple times the Cheb. polyn. and derivatives
!computed at (b,pk)
DO i=0,csi%order
    T1(i)=cos(i*acos(a1))
    T2(i)=cos(i*acos(a2))
END DO

val=0.

DO i1=0,csi%order
DO i2=0,csi%order
    val=val+coef(i1,i2)*T1(i1)*T2(i2)
END DO
END DO

END SUBROUTINE val0_chebyshev




SUBROUTINE val1_chebyshev(csi,b,pk,coef,val,d_wrt_b,d_wrt_pk)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Computes the value and the partial first order derivatives
!of the Chebyshev approximation at the point (b,pk)
!given the Chebyshev coefficients coef(0:order,0:order)
!___________________________________________________________ 
IMPLICIT NONE
TYPE (csi_type), intent(in)::csi
DOUBLE PRECISION, intent(in)::b,pk
DOUBLE PRECISION, dimension(0:csi%order,0:csi%order), intent(in)::coef
DOUBLE PRECISION, intent(out):: val,d_wrt_b,d_wrt_pk
INTEGER i,i1,i2
DOUBLE PRECISION a1,a2,d1,d2,c1,c2
DOUBLE PRECISION, dimension(0:csi%order)::T1,T2,T1x,T2x

IF ((b.GT.(csi%bmax)) .OR. (b.LT.(csi%bmin))) &
	STOP 'ERROR: interpolation bounds for b exceeded'

IF ((pk.GT.(csi%pkmax)) .OR. (pk.LT.(csi%pkmin))) &
	STOP 'ERROR: interpolation bounds for pk exceeded'

!The arguments of the polynomial and their derivative, and the 
!derivative coming from arccos
d1=2./(csi%bmax-csi%bmin)
d2=2./(csi%pkmax-csi%pkmin)
a1=(b-csi%bmin)*d1-1.
a2=(pk-csi%pkmin)*d2-1.
c1=1./sqrt(1.-a1**2)
c2=1./sqrt(1.-a2**2)

!Store once and use multiple times the Cheb. polyn. and derivatives
!computed at (b,pk)
DO i=0,csi%order
    T1(i)=cos(i*acos(a1))
    T2(i)=cos(i*acos(a2))
    T1x(i)=i*sin(i*acos(a1))
    T2x(i)=i*sin(i*acos(a2))
END DO

val=0.;d_wrt_b=0.;d_wrt_pk=0.

DO i1=0,csi%order
DO i2=0,csi%order
    val=val+coef(i1,i2)*T1(i1)*T2(i2)
    d_wrt_b=d_wrt_b+coef(i1,i2)*T1x(i1)*T2(i2)
    d_wrt_pk=d_wrt_pk+coef(i1,i2)*T1(i1)*T2x(i2)
END DO
END DO
d_wrt_b=d_wrt_b*d1*c1
d_wrt_pk=d_wrt_pk*d2*c2

END SUBROUTINE val1_chebyshev




SUBROUTINE val2_chebyshev(csi,b,pk,coef,val,E1,E2,E11,E12,E22)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Computes the value, first and second order partial derivatives
!of the Chebyshev approximation at the point (b,pk)
!given the Chebyshev coefficients coef(0:order,0:order)
!___________________________________________________________ 
IMPLICIT NONE
TYPE (csi_type), intent(in)::csi
DOUBLE PRECISION, intent(in)::b,pk
DOUBLE PRECISION, dimension(0:csi%order,0:csi%order), intent(in)::coef
DOUBLE PRECISION, intent(out):: val,E1,E2,E11,E12,E22
INTEGER i,i1,i2
DOUBLE PRECISION a1,a2,d1,d2,c1,c2
DOUBLE PRECISION, dimension(0:csi%order)::T1,T2,T1x,T2x,T1xx,T2xx

IF ((b.GT.(csi%bmax)) .OR. (b.LT.(csi%bmin))) &
	STOP 'ERROR: interpolation bounds for b exceeded'

IF ((pk.GT.(csi%pkmax)) .OR. (pk.LT.(csi%pkmin))) &
	STOP 'ERROR: interpolation bounds for pk exceeded'

!The arguments of the polynomial and their derivative, and the 
!derivative coming from arccos

d1=2./(csi%bmax-csi%bmin)
d2=2./(csi%pkmax-csi%pkmin)
a1=(b-csi%bmin)*d1-1.
a2=(pk-csi%pkmin)*d2-1.
c1=1./sqrt(1.-a1**2)
c2=1./sqrt(1.-a2**2)

!Store once and use multiple times the Cheb. polyn. and derivatives
!computed at (b,pk)
DO i=0,csi%order
    T1(i)=cos(i*acos(a1))
    T2(i)=cos(i*acos(a2))
    T1x(i)=i*sin(i*acos(a1))
    T2x(i)=i*sin(i*acos(a2))
    T1xx(i)=i*i*T1(i)
    T2xx(i)=i*i*T2(i)
END DO

val=0.;E1=0.;E2=0.;E11=0.;E12=0.;E22=0.

DO i1=0,csi%order
DO i2=0,csi%order
    val=val+coef(i1,i2)*T1(i1)*T2(i2)
    E1=E1+coef(i1,i2)*T1x(i1)*T2(i2)
    E2=E2+coef(i1,i2)*T1(i1)*T2x(i2)
    E11=E11+coef(i1,i2)*T1xx(i1)*T2(i2)
    E22=E22+coef(i1,i2)*T1(i1)*T2xx(i2)
    E12=E12+coef(i1,i2)*T1x(i1)*T2x(i2)
END DO
END DO

E11=(E1*a1*c1-E11)*d1**2*c1**2
E22=(E2*a2*c2-E22)*d2**2*c2**2
E1=E1*d1*c1
E2=E2*d2*c2
E12=E12*d1*d2*c1*c2

END SUBROUTINE val2_chebyshev




SUBROUTINE val0x_chebyshev(csi,b,pk,coef,val)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Computes Chebyshev approximation of function value
!at the point (b,pk) given the Chebyshev coefficients
!coef(0:order,0:order)
!Contrary to val0_chebyshev this function calculates
!extrapolated value if (b,pk) is outside the grid
!___________________________________________________________ 
IMPLICIT NONE
TYPE (csi_type), intent(in)::csi
DOUBLE PRECISION, intent(in)::b,pk
DOUBLE PRECISION, dimension(0:csi%order,0:csi%order), intent(in)::coef
DOUBLE PRECISION, intent(out):: val
INTEGER i,i1,i2
DOUBLE PRECISION a1,a2,d1,d2
DOUBLE PRECISION, dimension(0:csi%order)::T1,T2


!The arguments of the polynomial and their derivative, and the 
!derivative coming from arccos
d1=2./(csi%bmax-csi%bmin)
d2=2./(csi%pkmax-csi%pkmin)
a1=(b-csi%bmin)*d1-1.
a2=(pk-csi%pkmin)*d2-1.

!Store once and use multiple times the Cheb. polyn. and derivatives
!computed at (b,pk)
T1(0)=1.; T1(1)=a1
T2(0)=1.; T2(1)=a2
DO i=2,csi%order
    T1(i)=2.*a1*T1(i-1)-T1(i-2)
    T2(i)=2.*a2*T2(i-1)-T2(i-2)
END DO

val=0.

DO i1=0,csi%order
DO i2=0,csi%order
    val=val+coef(i1,i2)*T1(i1)*T2(i2)
END DO
END DO

END SUBROUTINE val0x_chebyshev
END MODULE chebyshev
