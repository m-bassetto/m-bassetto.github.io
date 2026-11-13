!****************************************************************************
!
! Purpose:
!  The main subroutines
!	 initcsi - compute initial guess for csi
!	 solvecsi - solve for csi by iterations
!	 solvess - solve for steady state given csi
!
! Last changed by Vadym Lepetyuk, October 11, 2004
! Copyright by Marco Bassetto, Florin Bidian, Vadym Lepetyuk
! This code can be freely distributed and modified for research purposes only, 
! provided this copyright notice is included in the modified code. 
! Proper credit should be given in all publications arising from
! modifications of this code; this should include a citation of 
! "Politics and Efficiency of Separating Capital and Ordinary Government Budgets"
! by Marco Bassetto with Thomas J. Sargent
!
!****************************************************************************


SUBROUTINE initgrid(csi)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Given number of nodes, and grid boundaries (bmin,bmax,pkmin,pkmax)
!allocate memory for and compute the Chebyshev interpolation nodes
!and the integral of the adjusted Chebyshev polynomials
!
!(IMPORTANT) The following fields of csi shoud be defined
!when this subroutine is called:
!csi%nodes, csi%bmin,csi%bmax,csi%pkmin,csi%pkmax
!_____________________________________________________________
!DEC$ REAL:8
USE modelsetup
USE chebyshev
IMPLICIT NONE
TYPE (csi_type), INTENT(inout)::csi
INTEGER i

csi%order=csi%nodes-1

ALLOCATE(csi%bnodes(csi%nodes))
ALLOCATE(csi%pknodes(csi%nodes))
ALLOCATE(csi%var(n_var,csi%nodes,csi%nodes))
ALLOCATE(csi%roots(csi%nodes))
ALLOCATE(csi%T(0:csi%order,1:csi%nodes))
ALLOCATE(csi%T2sum(0:csi%order))
ALLOCATE(csi%coef(n_var,0:csi%order,0:csi%order))

DO i=1,csi%nodes
    csi%roots(i)=-cos((2*i-1.)/(2*csi%nodes)*PIconst)  !roots solve  cos(nodes*acos(roots))=0
END DO

DO i=0,csi%order
    csi%T(i,:)=cos(i*acos(csi%roots))
END DO

DO i=0,csi%order
    csi%T2sum(i)=SUM(csi%T(i,:)**2)
END DO

ENTRY reinitgrid(csi)

csi%bnodes=(csi%roots+1)*(csi%bmax-csi%bmin)/2.+csi%bmin
csi%pknodes=(csi%roots+1)*(csi%pkmax-csi%pkmin)/2.+csi%pkmin

END SUBROUTINE initgrid




SUBROUTINE cleargrid(csi)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Release memory for all allocated arrays in csi
!_____________________________________________________________
!DEC$ REAL:8
USE modelsetup
TYPE (csi_type), intent(inout)::csi

DEALLOCATE(csi%bnodes)
DEALLOCATE(csi%pknodes)
DEALLOCATE(csi%var)
DEALLOCATE(csi%roots)
DEALLOCATE(csi%T)
DEALLOCATE(csi%T2sum)
DEALLOCATE(csi%coef)

END SUBROUTINE cleargrid




SUBROUTINE convertcsi(csi1,csi2)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Convert csi from one grid to the other.
! Extrapolation is used if necessary.
!_____________________________________________________________
!DEC$ REAL:8
USE modelsetup
USE chebyshev
IMPLICIT NONE
TYPE (csi_type), intent(in)::csi1
TYPE (csi_type), intent(inout)::csi2
INTEGER i,j1,j2

!Compute interpolation values (csi2%var) based on (csi1%coef)
DO j1=1,csi2%nodes
DO j2=1,csi2%nodes
	DO i=1,n_var
		CALL val0x_chebyshev(csi1,csi2%bnodes(j1),csi2%pknodes(j2),&
			csi1%coef(i,:,:),csi2%var(i,j1,j2))
	END DO
END DO
END DO

!Compute Chebyshev coefficients for all variables
!(this step is not really necessary)
!DO i=1,n_var
!	CALL coef_chebyshev(csi2,csi2%var(i,:,:),csi2%coef(i,:,:))
!END DO

END SUBROUTINE convertcsi




SUBROUTINE adjustboundaries(csi)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! (Suppose to be) sophisticated algorithm of adjusting
! boundaries of the grid (i.e. bmin,bmax,pkmin,pkmax)
!_____________________________________________________________
!DEC$ REAL:8
USE modelsetup
IMPLICIT NONE
TYPE (csi_type), intent(inout)::csi
DOUBLE PRECISION :: residbmin,residbmax,residpkmin,residpkmax

WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'adjustboundaries:Adjusting grid boundaries...'
WRITE(unit1,*) '************************************'
WRITE(unit1,*) '(bmin,bmax,pkmin,pkmax)'
WRITE(unit1,'(A12,4F10.3)') 'Old box', csi%bmin,csi%bmax,csi%pkmin,csi%pkmax
WRITE(unit1,'(A12,4F10.3)') '  csi  ', csi%var(i__b,1,csi%nodes),csi%var(i__b,csi%nodes,1),&
	csi%var(i__pk,csi%nodes,1),csi%var(i__pk,1,csi%nodes)

residbmin = (csi%var(i__b,1,csi%nodes)-csi%bmin)/(csi%bmax-csi%bmin)
residbmax = (csi%bmax-csi%var(i__b,csi%nodes,1))/(csi%bmax-csi%bmin)
residpkmin = (csi%var(i__pk,csi%nodes,1)-csi%pkmin)/(csi%pkmax-csi%pkmin)
residpkmax = (csi%pkmax-csi%var(i__pk,1,csi%nodes))/(csi%pkmax-csi%pkmin)
WRITE(unit1,'(A12,4F10.3)') ' resid ', residbmin,residbmax,residpkmin,residpkmax

!Box squeeze
 if (residbmin<adjtol) csi%pkmax=csi%pkmax-adj*(csi%pkmax-csi%pkmin)
!if (residbmax<adjtol) csi%pkmin=csi%pkmin+adj*(csi%pkmax-csi%pkmin)
!if (residpkmin<adjtol) csi%bmax=csi%bmax-adj*(csi%bmax-csi%bmin)
!if (residpkmax<adjtol) csi%bmin=csi%bmin+adj*(csi%bmax-csi%bmin)
WRITE(unit1,'(A12,4F10.3)') 'New box', csi%bmin,csi%bmax,csi%pkmin,csi%pkmax

END SUBROUTINE adjustboundaries




SUBROUTINE initcsi(csi)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Initial guess for csi values
!______________________________________________________________
!DEC$ REAL:8
USE modelsetup
IMPLICIT NONE
TYPE (csi_type), INTENT(inout)::csi
INTEGER :: j1,j2
INTEGER iostatus

SELECT CASE(999)

CASE(18)  !reading csi from file
	CALL readcsi18('csi18.txt',csi)

CASE(999)	 !sophisticated guess
	csi%var=0.
	DO j1=1,csi%nodes	!analytical solution for log preferences
	DO j2=1,csi%nodes	!for the case with no public sector (i.e. G=pk=0)
		csi%var(i__gam,j1,j2)=-(1.-deltapk)/(1.+n)*csi%pknodes(j2)
		csi%var(i__tau,j1,j2)=(1.+alpha*p)/(1.+theta+n)*csi%bnodes(j1)&
			+(1.-x)*qpk*csi%var(i__gam,j1,j2)
		csi%var(i__b,j1,j2)=(1.-alpha)/(1.+n)*csi%bnodes(j1)&
			+(1.+theta+n)/(1.+n)*x/p*qpk*csi%var(i__gam,j1,j2)
		csi%var(i__cy,j1,j2)=(&
			w-(1+n+theta*x)/(1.+n)*qpk*csi%var(i__gam,j1,j2)&
			+( qk/(r+(1-deltak)*qk)*csi%var(i__b,j1,j2)-csi%bnodes(j1)/(1.+n) )&
				*( 1.+p-theta*(1.+alpha*p)/(1.+theta+n) )&
			) / (2.+beta*theta)
		csi%var(i__l,j1,j2)=1.-csi%var(i__cy,j1,j2)/w
		csi%var(i__k,j1,j2)=( w*csi%var(i__l,j1,j2)-csi%var(i__cy,j1,j2)&
			-((p+1.)*csi%bnodes(j1)-theta*csi%var(i__tau,j1,j2))/(1.+n)&
			-(1.+theta+n)/(1.+n)*qpk*csi%var(i__gam,j1,j2) )/qk
	END DO
	END DO
	csi%var(i__G,:,:)=0.2
	csi%var(i__pk,:,:)=(csi%pkmax+csi%pkmin)/2.

CASE DEFAULT	 !simple guess
	csi%var=0.
	csi%var(i__cy,:,:)=1.
	DO j1=1,csi%nodes
	DO j2=1,csi%nodes
		!implicit assumption: G=pk=1.
		csi%var(i__tau,j1,j2)=(1.+alpha*p)/(1.+theta+n)*csi%bnodes(j1)+&
			(1.-y)*qG+(1.-x)*qpk*(1.-(1.-deltapk)/(1.+n)*csi%pknodes(j2))
	END DO
	END DO
	csi%var(i__l,:,:)=.5
	csi%var(i__k,:,:)=1.
	csi%var(i__b,:,:)=(csi%bmax+csi%bmin)/2.
	csi%var(i__G,:,:)=1.
	csi%var(i__pk,:,:)=(csi%pkmax+csi%pkmin)/2.

END SELECT
END SUBROUTINE initcsi




SUBROUTINE testcsi(csi)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Initial tests (optional)
!______________________________________________________________
!DEC$ REAL:8
USE modelsetup !Set up the model and global variables and routines
USE chebyshev !Chebyshev approximation subroutines
USE numerical_libraries !Defines IMSL's routine interfaces
IMPLICIT NONE
TYPE (csi_type), INTENT(inout), TARGET::csi
EXTERNAL fcn,lsjac
INTEGER, PARAMETER:: ldfjac=n_eqn
DOUBLE PRECISION, DIMENSION(n_eqn):: fvec
DOUBLE PRECISION, DIMENSION(n_var):: xlb, xub
DOUBLE PRECISION, DIMENSION(n_eqn,n_var) :: fjac,fjac2
DOUBLE PRECISION :: rparam(7), rp(7) ! IMSL BCLSJ parameters
INTEGER :: ibtype, iparam(6), ip(6) ! IMSL BCLSJ parameters
DOUBLE PRECISION, DIMENSION(n_eqn):: fscale ! IMSL BCLSJ parameters
DOUBLE PRECISION, DIMENSION(n_var):: xscale ! IMSL BCLSJ parameters
DOUBLE PRECISION, DIMENSION(n_var):: xguess,xx
INTEGER, DIMENSION(n_eqn,n_var) :: infojac
INTEGER, DIMENSION(2):: maxlocjac


! IMSL BCLSJ parameters
ibtype=0   !User will supply all the bounds.
fscale=1.0 !Vector of the diagonal scaling matrix for the functions
xscale=1.0 !Vector of the diagonal scaling matrix for the variables

xlb=-1.0E30
xlb(i__cy)=.0001
xlb(i__l)=.0001
xlb(i__b)=csi%bmin+.0001 ! b: this bound due to interpolation
xlb(i__G)=.0001
xlb(i__pk)=csi%pkmin+.0001 ! pk: this bound due to interpolation

xub=+1.0E30
xub(i__l)=0.9999
xub(i__b)=csi%bmax-.0001 ! b : this limit due to interpolation
xub(i__pk)=csi%pkmax-.0001 ! pk: this limit due to interpolation

!IMSL:U4LSF Provide default parameters for the nonlinear least squares problems.
CALL DU4LSF(iparam,rparam)
iparam(1)=1 !Initialization flag (1 for user-specified parameters)
iparam(3)=100000 !Maximum number of iterations
iparam(4)=100000 !Maximum number of function evaluations
iparam(5)=100000 !Maximum number of Jacobian evaluations
rparam(3)=0.001*rparam(3) !Relative function tolerance


WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'main:Initial Tests...'
WRITE(unit1,*) '************************************'
!Current csi vector
o_csi=>csi
!Current proxies for tau,G,pk
CALL coef_chebyshev(csi,csi%var(i__tau,:,:),csi%coef(i__tau,:,:))
CALL coef_chebyshev(csi,csi%var(i__G,:,:),csi%coef(i__G,:,:))
CALL coef_chebyshev(csi,csi%var(i__pk,:,:),csi%coef(i__pk,:,:))
!Current state of the economy
o_b=csi%bnodes(4);o_pk=csi%pknodes(4)
!Initial guess
xguess=csi%var(:,4,4)

WRITE(unit1,*) 'testcsi:Checking fcn...'
CALL writexx(unit1,'testcsi:xguess',n_var,xguess)
CALL fcn(n_eqn,n_eqn,xguess,fvec)
CALL writeff(unit1,'testcsi:f(xguess)',n_eqn,fvec)

!   IMSL:FDJAC Approximate the Jacobian using forward differences 
WRITE(unit1,*) 'testcsi:Comparing finite-difference and user-supplied Jacobians...'
CALL DFDJAC(fcn,n_eqn,n_eqn,xguess,xscale,fvec,0.0,fjac,n_eqn)
CALL lsjac(n_eqn,n_eqn,xguess,fjac2,n_eqn)
WRITE(unit1,*) 'MaxAbs Dirrefence', maxval(abs(fjac2-fjac))
maxlocjac=maxloc(abs(fjac2-fjac))
WRITE(unit1,*) 'MaxAbs Location', maxlocjac
WRITE(unit1,*) 'MaxAbs Finite Differ Jacobian', fjac(maxlocjac(1),maxlocjac(2))
WRITE(unit1,*) 'MaxAbs User-supplied Jacobian', fjac2(maxlocjac(1),maxlocjac(2))

!   IMSL:BCLSF Solve a nonlinear least squares problem subject to bounds
!   on the variables using a modified Levenberg-Marquardt algorithm
!   and a finite-difference Jacobian.
WRITE(unit1,*) 'testcsi:Calling BCLSF (finite-difference Jacobian)...'
ip=iparam;rp=rparam
CALL DBCLSF(fcn,n_eqn,n_var,xguess,ibtype,xlb,xub,xscale,&
	  &fscale,ip,rp,xx,fvec,fjac,ldfjac)
WRITE(unit1,*) 'testcsi:Getting results from BCLSF...'
CALL writepx(unit1,'testcsi:RESULT',ip)
CALL writexx(unit1,'testcsi:xx',n_var,xx)
CALL writeff(unit1,'testcsi:f(xx)',n_eqn,fvec)

!   IMSL:BCLSJ Solve a nonlinear least squares problem subject to bounds
!   on the variables using a modified Levenberg-Marquardt algorithm
!   and a user-supplied Jacobian.
WRITE(unit1,*) 'testcsi:Calling BCLSJ (user-supplied Jacobian)...'
ip=iparam;rp=rparam
CALL DBCLSJ(fcn,lsjac,n_eqn,n_eqn,xguess,ibtype,xlb,xub,xscale,&
	  &fscale,ip,rp,xx,fvec,fjac,ldfjac)
WRITE(unit1,*) 'testcsi:Getting results from BCLSJ...'
CALL writepx(unit1,'testcsi:RESULT',ip)
CALL writexx(unit1,'testcsi:xx',n_var,xx)
CALL writeff(unit1,'testcsi:f(xx)',n_eqn,fvec)

END SUBROUTINE testcsi




SUBROUTINE solvecsi(csi)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Solve for csi
!______________________________________________________________
!DEC$ REAL:8
USE modelsetup !Set up the model and global variables and routines
USE chebyshev !Chebyshev approximation subroutines
USE numerical_libraries !Defines IMSL's routine interfaces
IMPLICIT NONE
TYPE (csi_type), INTENT(inout), TARGET::csi
EXTERNAL fcn,lsjac
INTEGER, PARAMETER:: ldfjac=n_eqn
DOUBLE PRECISION, DIMENSION(n_eqn):: fvec
DOUBLE PRECISION, DIMENSION(n_var):: xlb, xub
DOUBLE PRECISION, DIMENSION(n_eqn,n_var) :: fjac
DOUBLE PRECISION :: rparam(7), rp(7) ! IMSL BCLSJ parameters
INTEGER :: ibtype, iparam(6), ip(6) ! IMSL BCLSJ parameters
DOUBLE PRECISION, DIMENSION(n_eqn):: fscale ! IMSL BCLSJ parameters
DOUBLE PRECISION, DIMENSION(n_var):: xscale ! IMSL BCLSJ parameters
DOUBLE PRECISION, DIMENSION(n_var):: xguess,xx
INTEGER :: it,i,j1,j2
DOUBLE PRECISION, PARAMETER :: tol=0.000001 ! Convergency tolerance
INTEGER, PARAMETER :: maxit=1000 ! Maximum number of iterations
DOUBLE PRECISION, DIMENSION(n_var,csi%nodes,csi%nodes)::newcsivar
DOUBLE PRECISION loopdiff
INTEGER iostatus


WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'solvecsi:Solving for csi...'
WRITE(unit1,*) '************************************'
CALL writesetup(unit1,'solvecsi:Parameters',csi)

! IMSL BCLSJ parameters
ibtype=0   !User will supply all the bounds.
fscale=1.0 !Vector of the diagonal scaling matrix for the functions
xscale=1.0 !Vector of the diagonal scaling matrix for the variables

xlb=-1.0E30
xlb(i__cy)=.0001
xlb(i__l)=.0001
xlb(i__b)=csi%bmin+.0001 ! b: this bound due to interpolation
xlb(i__G)=.0001
xlb(i__pk)=csi%pkmin+.0001 ! pk: this bound due to interpolation

xub=+1.0E30
xub(i__l)=0.9999
xub(i__b)=csi%bmax-.0001 ! b : this limit due to interpolation
xub(i__pk)=csi%pkmax-.0001 ! pk: this limit due to interpolation

!IMSL:U4LSF Provide default parameters for the nonlinear least squares problems.
CALL DU4LSF(iparam,rparam)
iparam(1)=1 !Initialization flag (1 for user-specified parameters)
iparam(3)=100000 !Maximum number of iterations
iparam(4)=100000 !Maximum number of function evaluations
iparam(5)=100000 !Maximum number of Jacobian evaluations
rparam(3)=0.001*rparam(3) !Relative function tolerance

CALL writepp(unit1,'solvecsi:BCLSJ',iparam,rparam)


!WRITE(unit1,*) '************************************'
!WRITE(unit1,*) 'solvecsi:Main Loop...'
!WRITE(unit1,*) '************************************'
it=0
DO WHILE (.TRUE.) !the loop is ended by EXIT statement
	it=it+1
	WRITE(unit1,*) 'solvecsi:Iteration ', it
	WRITE(unit1,*) '(it, j1,j2, b,pk, ssr)'

	o_csi=>csi
	CALL coef_chebyshev(csi,csi%var(i__tau,:,:),csi%coef(i__tau,:,:))
	CALL coef_chebyshev(csi,csi%var(i__G,:,:),csi%coef(i__G,:,:))
	CALL coef_chebyshev(csi,csi%var(i__pk,:,:),csi%coef(i__pk,:,:))

	DO j1=1,csi%nodes
	DO j2=1,csi%nodes

		o_b=csi%bnodes(j1); o_pk=csi%pknodes(j2)
		xguess=csi%var(:,j1,j2)

		ip=iparam;rp=rparam
		CALL DBCLSJ(fcn,lsjac,n_eqn,n_var,xguess,ibtype,xlb,xub,xscale,&
	!   CALL DBCLSF(fcn,n_eqn,n_var,xguess,ibtype,xlb,xub,xscale,&
			&fscale,ip,rp,xx,fvec,fjac,ldfjac)

		WRITE(unit1,'(3I4,2F12.5,E14.4)') it,j1,j2, o_b,o_pk, SUM(fvec**2)
		CALL writeboundhit(unit1,n_var,xx,xlb,xub)
		
		newcsivar(:,j1,j2)=xx

	END DO
	END DO

	loopdiff=maxval(abs(csi%var-newcsivar))
	WRITE(unit1,*) 'solvecsi:Iteration Difference ', loopdiff

	csi%var=newcsivar

	IF(loopdiff<tol.OR.it>=maxit) EXIT !from DO WHILE loop

END DO !WHILE
WRITE(unit1,*) 'solvecsi:Number of iterations ', it
IF (it>=maxit) WRITE(unit1,*) 'solvecsi: WARNING: Maximum number of iterations exceeded'

!Compute Chebyshev coefficients for all variables
DO i=1,n_var
	CALL coef_chebyshev(csi,csi%var(i,:,:),csi%coef(i,:,:))
END DO


!Writing csi into a text file
CALL writecsi(csi)

!Writing csi into Array Visualizer's files
CALL writeagl(csi)

END SUBROUTINE solvecsi




SUBROUTINE solvess(csi,xss)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Solve for steady state given csi
!______________________________________________________________
!DEC$ REAL:8
USE modelsetup
USE chebyshev
IMPLICIT NONE
TYPE (csi_type), INTENT(in), TARGET::csi
DOUBLE PRECISION, DIMENSION(n_var), INTENT(out):: xss
DOUBLE PRECISION, DIMENSION(n_var):: dxssb,dxsspk
DOUBLE PRECISION ::loopdiff
INTEGER ::i,it
DOUBLE PRECISION, PARAMETER :: tol=0.000001 ! Convergency tolerance
INTEGER, PARAMETER :: maxit=1000 ! Maximum number of iterations
DOUBLE PRECISION ::b,pk,b1,pk1

WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'solvess:Looking for steady state ...'
WRITE(unit1,*) '************************************'

b=csi%bnodes((csi%nodes+1)/2)
pk=csi%pknodes((csi%nodes+1)/2)
!WRITE(unit1,'(A12,2F12.5)') 'b, pk', b, pk

it=0
DO WHILE (.TRUE.) !the loop is ended by EXIT statement
	it=it+1
	!WRITE(unit1,*) 'solvess:Iteration ', it

	o_csi=>csi
	CALL val0_chebyshev(csi,b,pk,csi%coef(i__b,:,:),b1)
	CALL val0_chebyshev(csi,b,pk,csi%coef(i__pk,:,:),pk1)

	loopdiff=max(abs(b-b1),abs(pk-pk1))
	b=b1;pk=pk1;
!	WRITE(unit1,'(A12,2F12.5)') 'b, pk', b, pk

	IF(loopdiff<tol.OR.it>=maxit) EXIT !from DO WHILE loop

END DO !WHILE
WRITE(unit1,*) 'solvess:Number of iterations ', it
IF (it>=maxit) WRITE(unit1,*) 'solvess: WARNING: Maximum number of iterations exceeded'

!Compute all variables at the steady state
DO i=1,n_var
	CALL val1_chebyshev(csi,b,pk,csi%coef(i,:,:),xss(i),dxssb(i),dxsspk(i))
END DO

CALL writexxs(unit1,'solvess: Steady State',n_var,xss)
WRITE(unit1,'(A16,2F12.5)') 'dv of future b ', dxssb(i__b), dxsspk(i__b)
WRITE(unit1,'(A16,2F12.5)') 'dv of future pk', dxssb(i__pk), dxsspk(i__pk)
WRITE(unit1,'(A16,2F12.5)') 'dv of future G ', dxssb(i__G), dxsspk(i__G)
END SUBROUTINE solvess




SUBROUTINE compensation(csi,xx0,mc)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Solve for a compensation for regime change
! (only valid for log utility)
!______________________________________________________________
!DEC$ REAL:8
USE modelsetup
USE chebyshev
IMPLICIT NONE
TYPE (csi_type), INTENT(in), TARGET::csi
DOUBLE PRECISION, DIMENSION(n_var), INTENT(in):: xx0
DOUBLE PRECISION, INTENT(out):: mc
DOUBLE PRECISION, DIMENSION(n_var) :: xx, xxn
DOUBLE PRECISION :: u, vxx, uxx0, m1,cotemp,transwedge,pvt
DOUBLE PRECISION ::loopdiff
INTEGER ::i,it
DOUBLE PRECISION, PARAMETER :: tol=0.000001 ! Convergency tolerance
INTEGER, PARAMETER :: maxit=1000 ! Maximum number of iterations

WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'compensation:Computing regime change compensation...'
WRITE(unit1,*) '************************************'

u=utility(xx0,xx0)
WRITE(unit1,*) '(b, pk, u, m1, mc, pvt, transwedge)'

pvt=xx0(i__tau)*(1.+theta/(1.+r-deltak)) !present value of taxes
uxx0=(2.+tbeta)*log(w-pvt)+&
	(1.+tbeta)*(phi*log(xx0(i__G))+eta*log(xx0(i__pk))) !steady state utility
	!excluding the constant:
	!-log(w)+tbeta*log(beta*(1.+r-deltak))-(2.+tbeta)*log(2.+tbeta)
WRITE(unit1,'(6F12.5)') xx0(i__b),xx0(i__pk),u,0.0,0.0,pvt

xx=xx0
cotemp=((r+(1.-deltak))*qk*xx0(i__k)+(p+1.)*xx0(i__b))/theta-xx0(i__tau)
mc=0.
it=0
DO WHILE (.TRUE.) !the loop is ended by EXIT statement
	it=it+1
	!WRITE(unit1,*) 'compensation:Iteration ', it

	DO i=1,n_var
		CALL val0_chebyshev(csi,xx(i__b),xx(i__pk),csi%coef(i,:,:),xxn(i))
	END DO

	u=utility(xx,xxn)
	vxx=phi*log(xx(i__G))+eta*log(xx(i__pk))+&
		tbeta*(phi*log(xxn(i__G))+eta*log(xxn(i__pk)))
	pvt=xx(i__tau)+theta*xxn(i__tau)/(1.+r-deltak) !present value of taxes
	m1=exp((uxx0-vxx)/(2.+tbeta))-w+pvt !one-period compensation
	IF ((w-pvt+m1)/(2.+tbeta)>w) WRITE(unit1,*) 'WARNING: l is negative'
	mc=mc+m1*((1.+n)/(1.+r-deltak))**(it-1) !cumulative discounted compensation
	transwedge=v_G(xx(i__G))*(xx(i__cy)+theta*cotemp/(1.+n))/(1.+theta/(1.+n))

	loopdiff=max(abs(xx(i__b)-xxn(i__b)),abs(xx(i__pk)-xxn(i__pk)))
	cotemp=((r+(1.-deltak))*qk*xx(i__k)+(p+1.)*xx(i__b))/theta-xxn(i__tau)
	xx=xxn
	WRITE(unit1,'(7F12.5)') xx(i__b),xx(i__pk),u,m1,mc,pvt,transwedge

	IF(loopdiff<tol.OR.it>=maxit) EXIT !from DO WHILE loop

END DO !WHILE
!WRITE(unit1,*) 'compensation:Number of iterations ', it
IF (it>=maxit) WRITE(unit1,*) 'compensation: WARNING: Maximum number of iterations exceeded'

!Long-run compensation
u=utility(xx,xx)
vxx=(1.+tbeta)*(phi*log(xx(i__G))+eta*log(xx(i__pk)))
pvt=xx(i__tau)*(1.+theta/(1.+r-deltak)) !present value of taxes
m1=exp((uxx0-vxx)/(2.+tbeta))-w+pvt !one-period compensation
mc=mc+m1*((1.+n)/(1.+r-deltak))**(it)/(1.-(1.+n)/(1.+r-deltak)) !cumulative discounted compensation
cotemp=((r+(1.-deltak))*qk*xx(i__k)+(p+1.)*xx(i__b))/theta-xxn(i__tau)
transwedge=v_G(xx(i__G))*(xx(i__cy)+theta*cotemp/(1.+n))/(1.+theta/(1.+n))
WRITE(unit1,'(7F12.5)') xx(i__b),xx(i__pk),u,m1,mc,pvt,transwedge

WRITE(unit1,*) 'compensation:Compensation', mc

END SUBROUTINE compensation




SUBROUTINE compensation1(csi,xx0,mc)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Solve for a compensation for regime change
! (only valid for log utility)
!______________________________________________________________
!DEC$ REAL:8
USE modelsetup
USE chebyshev
IMPLICIT NONE
TYPE (csi_type), INTENT(in), TARGET::csi
DOUBLE PRECISION, DIMENSION(n_var), INTENT(in):: xx0
DOUBLE PRECISION, INTENT(out):: mc
DOUBLE PRECISION, DIMENSION(n_var) :: xx, xxn
DOUBLE PRECISION :: u, uxx, vxx0, pvt, pvt0, m1
DOUBLE PRECISION ::loopdiff
INTEGER ::i,it
DOUBLE PRECISION, PARAMETER :: tol=0.000001 ! Convergency tolerance
INTEGER, PARAMETER :: maxit=1000 ! Maximum number of iterations

WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'compensation:Computing regime change compensation...'
WRITE(unit1,*) '************************************'

u=utility(xx0,xx0)
WRITE(unit1,*) '(b, pk, u, m1, mc, pvt)'

vxx0=(1.+tbeta)*(phi*log(xx0(i__G))+eta*log(xx0(i__pk)))
pvt0=xx0(i__tau)*(1.+theta/(1.+r-deltak)) !present value of taxes
WRITE(unit1,'(6F12.5)') xx0(i__b),xx0(i__pk),u,0.0,0.0,pvt0

xx=xx0
mc=0.
it=0
DO WHILE (.TRUE.) !the loop is ended by EXIT statement
	it=it+1
	!WRITE(unit1,*) 'compensation:Iteration ', it

	DO i=1,n_var
		CALL val0_chebyshev(csi,xx(i__b),xx(i__pk),csi%coef(i,:,:),xxn(i))
	END DO

	u=utility(xx,xxn)
	pvt=xx(i__tau)+theta*xxn(i__tau)/(1.+r-deltak) !present value of taxes
	uxx=(2.+tbeta)*log(w-pvt)+&
		phi*log(xx(i__G))+eta*log(xx(i__pk))+&
		tbeta*(phi*log(xxn(i__G))+eta*log(xxn(i__pk)))
	m1=w-pvt0-exp((uxx-vxx0)/(2.+tbeta)) !one-period compensation
	IF ((w-pvt-m1)/(2.+tbeta)>w) WRITE(unit1,*) 'WARNING: l is negative'
	mc=mc+m1*((1.+n)/(1.+r-deltak))**(it-1) !cumulative discounted compensation
	
	loopdiff=max(abs(xx(i__b)-xxn(i__b)),abs(xx(i__pk)-xxn(i__pk)))
	xx=xxn
	WRITE(unit1,'(6F12.5)') xx(i__b),xx(i__pk),u,m1,mc,pvt

	IF(loopdiff<tol.OR.it>=maxit) EXIT !from DO WHILE loop

END DO !WHILE
!WRITE(unit1,*) 'compensation:Number of iterations ', it
IF (it>=maxit) WRITE(unit1,*) 'compensation: WARNING: Maximum number of iterations exceeded'

!Long-run compensation
u=utility(xx,xx)
pvt=xx(i__tau)*(1.+theta/(1.+r-deltak)) !present value of taxes
uxx=(2.+tbeta)*log(w-pvt)+&
	(1.+tbeta)*(phi*log(xx(i__G))+eta*log(xx(i__pk)))
m1=w-pvt0-exp((uxx-vxx0)/(2.+tbeta)) !one-period compensation
mc=mc+m1*((1.+n)/(1.+r-deltak))**(it)/(1.-(1.+n)/(1.+r-deltak)) !cumulative discounted compensation
WRITE(unit1,'(6F12.5)') xx(i__b),xx(i__pk),u,m1,mc,pvt

WRITE(unit1,*) 'compensation:Compensation', mc

END SUBROUTINE compensation1




SUBROUTINE efficiency(xss,wedgeG,wedgepk)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Solve for asymptotic efficiency
!______________________________________________________________
!DEC$ REAL:8
USE modelsetup
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(n_var), INTENT(in):: xss
DOUBLE PRECISION, INTENT(out):: wedgeG,wedgepk
DOUBLE PRECISION:: co

WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'efficiency:Computing wedges...'
WRITE(unit1,*) '************************************'

co=((r+(1.-deltak))*qk*xss(i__k)+(p+1.)*xss(i__b))/theta-xss(i__tau)

wedgeG=v_G(xss(i__G))*(1.+n)/(1.+theta+n)*&
	(1./u_y_c(xss(i__cy),xss(i__l))+theta/((1.+n)*u_o_c(co)))
wedgepk=v_pk(xss(i__pk))*(1.+n)/(1.+theta+n)*&
	(1./u_y_c(xss(i__cy),xss(i__l))+theta/((1.+n)*u_o_c(co)))&
	/(1.-(1.-deltapk)/(1.+r-deltapk))
WRITE(unit1,*) 'efficiency:wedgeG', wedgeG
WRITE(unit1,*) 'efficiency:wedgepk', wedgepk

END SUBROUTINE efficiency
