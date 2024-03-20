!****************************************************************************
!
! Purpose:
!  A few subroutines solely to generate "nice" output
!
! Last changed by Vadym Lepetyuk, March 27, 2005
! Copyright by Marco Bassetto, Florin Bidian, Vadym Lepetyuk
! This code can be freely distributed and modified for research purposes only, 
! provided this copyright notice is included in the modified code. 
! Proper credit should be given in all publications arising from
! modifications of this code; this should include a citation of 
! "Politics and Efficiency of Separating Capital and Ordinary Government Budgets"
! by Marco Bassetto with Thomas J. Sargent
!
!****************************************************************************

SUBROUTINE writexx(unit,title,nn,xx)
!DEC$ REAL:8
USE modelsetup
IMPLICIT NONE
INTEGER, INTENT(in):: unit
CHARACTER (LEN=*), INTENT(in):: title
INTEGER, INTENT(in):: nn
DOUBLE PRECISION, DIMENSION(nn), INTENT(in):: xx
DOUBLE PRECISION cy,l,gam,k,b,tau,G,pk,&
	cy_G,cy_pk,l_G,l_pk,k_G,k_pk,b_G,b_pk,&
	tau_G,tau_pk,w_G,w_pk,r_G,r_pk,p_G,p_pk

cy=xx(i__cy);l=xx(i__l);gam=xx(i__gam)
k=xx(i__k);b=xx(i__b);tau=xx(i__tau)
G=xx(i__G);pk=xx(i__pk);cy_G=xx(i__cy_G);cy_pk=xx(i__cy_pk)
l_G=xx(i__l_G);l_pk=xx(i__l_pk);k_G=xx(i__k_G);k_pk=xx(i__k_pk)
b_G=xx(i__b_G);b_pk=xx(i__b_pk);tau_G=xx(i__tau_G);tau_pk=xx(i__tau_pk)

WRITE(unit,*) '---------------- ',title,' ----------------'
WRITE(unit,'(6(A5,3F12.5/),(A5,F12.5/),(A5,F12.5))') &
	'cy',   cy,  cy_G,  cy_pk, &  ! Consumption when young
	'l',	l,   l_G,   l_pk, &   ! Labor when young
	'gam',  gam, 0.0,   1.0, &    ! Public Investment
	'k',	k,   k_G,   k_pk, &   ! Private Capital
	'b',	b,   b_G,   b_pk, &   ! Bonds
	'tau',  tau, tau_G, tau_pk, & ! Lump-sum Taxes
	'G',	G, &				  ! Nondurable Public Good
	'pk',   pk					  ! Durable Public Good

END SUBROUTINE writexx



SUBROUTINE writexxs(unit,title,nn,xss)
!DEC$ REAL:8
USE modelsetup
IMPLICIT NONE
INTEGER, INTENT(in):: unit
CHARACTER (LEN=*), INTENT(in):: title
INTEGER, INTENT(in):: nn
DOUBLE PRECISION, DIMENSION(nn), INTENT(in):: xss
DOUBLE PRECISION co

CALL writexx(unit,title,nn,xss)

co=((r+(1.-deltak))*qk*xss(i__k)+(p+1.)*xss(i__b))/theta-xss(i__tau)
WRITE(unit,'(A5,F12.5)') 'co',co  ! Steady state co

END SUBROUTINE writexxs



SUBROUTINE writeff(unit,title,nn,ff)
IMPLICIT NONE
INTEGER, INTENT(in):: unit
CHARACTER (LEN=*), INTENT(in):: title
INTEGER, INTENT(in):: nn
DOUBLE PRECISION, DIMENSION(nn), INTENT(in):: ff
DOUBLE PRECISION :: ssr
INTEGER :: i, i1, i2, i3

WRITE(unit,*) '---------------- ',title,' ----------------'
DO i=1,nn/3
	i1=i; i2=i+nn/3; i3=i+2*nn/3
	WRITE(unit,'(3("f(",I2,")",F12.5,"	"))') &
		i1,ff(i1),i2,ff(i2),i3,ff(i3)
END DO

ssr=SUM(ff**2)
WRITE(unit,*) '>>>>>>>>> ssr = ', ssr

END SUBROUTINE writeff



SUBROUTINE writepp(unit,title,iparam,rparam)
IMPLICIT NONE
INTEGER, INTENT(in):: unit
CHARACTER (LEN=*), INTENT(in):: title
INTEGER, INTENT(in):: iparam(6)
DOUBLE PRECISION, INTENT(in) :: rparam(7)

WRITE(unit,*) '---------------- ',title,':iparam ----------------'
WRITE(unit,'(A44,I12)') &
	'Initialization flag', iparam(1), &
	'Number of good digits in the function', iparam(2), &
	'Maximum number of iterations', iparam(3), &
	'Maximum number of function evaluations', iparam(4), &
	'Maximum number of Jacobian evaluations', iparam(5), &
	'Internal variable scaling flag', iparam(6)
WRITE(unit,*) '---------------- ',title,':rparam ----------------'
WRITE(unit,'(A44,G12.4)') &
	'Scaled gradient tolerance', rparam(1), &
	'Scaled step tolerance', rparam(2), &
	'Relative function tolerance', rparam(3), &
	'Absolute function tolerance', rparam(4), &
	'False convergence tolerance', rparam(5), &
	'Maximum allowable step size', rparam(6), &
	'Size of initial trust region radius', rparam(7)
END SUBROUTINE writepp



SUBROUTINE writepx(unit,title,iparam)
IMPLICIT NONE
INTEGER, INTENT(in):: unit
CHARACTER (LEN=*), INTENT(in):: title
INTEGER, INTENT(in):: iparam(6)

WRITE(unit,*) '---------------- ',title,' ----------------'
WRITE(unit,'(A34,I12)') &
	'Number of iterations', iparam(3), &
	'Number of function evaluations', iparam(4), &
	'Number of Jacobian evaluations', iparam(5)
END SUBROUTINE writepx



SUBROUTINE writesetup(unit,title,csi)
USE modelsetup
IMPLICIT NONE
INTEGER, INTENT(in):: unit
CHARACTER (LEN=*), INTENT(in):: title
TYPE (csi_type), INTENT(in)::csi

WRITE(unit,*) '---------------- ',title,' ----------------'
WRITE(unit,'(A44,F12.3)') 'Discount factor (beta)', beta
WRITE(unit,'(A44,F12.3)') 'Probability of surviving (theta)', theta
WRITE(unit,'(A44,F12.3)') 'Population growth (n)', n
WRITE(unit,'(A44,F12.3)') 'Debt repayment rate (alpha)', alpha
WRITE(unit,'(A44,F12.3)') 'Public investment paid by bonds (x)', x
WRITE(unit,'(A44,F12.3)') 'Public consumption paid by bonds (y)', y

WRITE(unit,*) 'Technology and shocks:'
WRITE(unit,'(A44,F12.3)') 'Depreciation rate for k', deltak
WRITE(unit,'(A44,F12.3)') 'Depreciation rate for pk', deltapk
WRITE(unit,'(A44,F12.3)') 'Relative price of capital', qk
WRITE(unit,'(A44,F12.3)') 'Relative price of public consumption', qG
WRITE(unit,'(A44,F12.3)') 'Relative price of public capital', qpk

WRITE(unit,*) 'Prices:'
WRITE(unit,'(A44,F12.3)') 'Return on capital (r)', r
WRITE(unit,'(A44,F12.3)') 'Return on labor (w)', w
WRITE(unit,'(A44,F12.3)') 'Price of console (p)', p

WRITE(unit,*) 'Preferences:'
!WRITE(unit,'(A44,F12.3)') 'sigma', sigma
!WRITE(unit,'(A44,F12.3)') 'miu', miu
WRITE(unit,'(A44,F12.3)') 'phi', phi
WRITE(unit,'(A44,F12.3)') 'eta', eta

WRITE(unit,*) 'Chebyshev interpolation:'
WRITE(unit,'(A44,I12)') 'Number of nodes', csi%nodes
WRITE(unit,'(A44,I12)') 'Order of polynomials', csi%order
WRITE(unit,'(A44,F12.3)') 'bmin', csi%bmin
WRITE(unit,'(A44,F12.3)') 'bmax', csi%bmax
WRITE(unit,'(A44,F12.3)') 'pkmin', csi%pkmin
WRITE(unit,'(A44,F12.3)') 'pkmax', csi%pkmax

END SUBROUTINE writesetup



SUBROUTINE writeboundhit(unit,nn,xx,xlb,xub)
!DEC$ REAL:8
USE modelsetup
IMPLICIT NONE
INTEGER, INTENT(in):: unit
INTEGER, INTENT(in):: nn
DOUBLE PRECISION, DIMENSION(nn), INTENT(in):: xx,xlb,xub
CHARACTER (LEN=80) :: str

str=''
IF (ABS((xx(i__cy)- xlb(i__cy))/ xlb(i__cy)) <.01) str=TRIM(str)//'cy- '
IF (ABS((xx(i__l)-  xlb(i__l))/  xlb(i__l))  <.01) str=TRIM(str)//'l- '
IF (ABS((xx(i__gam)-xlb(i__gam))/xlb(i__gam))<.01) str=TRIM(str)//'gam- '
IF (ABS((xx(i__k)-  xlb(i__k))/  xlb(i__k))  <.01) str=TRIM(str)//'k- '
IF (ABS((xx(i__b)-  xlb(i__b))/  xlb(i__b))  <.01) str=TRIM(str)//'b- '
IF (ABS((xx(i__tau)-xlb(i__tau))/xlb(i__tau))<.01) str=TRIM(str)//'tau- '
IF (ABS((xx(i__G)-  xlb(i__G))/  xlb(i__G))  <.01) str=TRIM(str)//'G- '
IF (ABS((xx(i__pk)- xlb(i__pk))/ xlb(i__pk)) <.01) str=TRIM(str)//'pk- '

IF (ABS((xx(i__cy)- xub(i__cy))/ xub(i__cy)) <.01) str=TRIM(str)//'cy+ '
IF (ABS((xx(i__l)-  xub(i__l))/  xub(i__l))  <.01) str=TRIM(str)//'l+ '
IF (ABS((xx(i__gam)-xub(i__gam))/xub(i__gam))<.01) str=TRIM(str)//'gam+ '
IF (ABS((xx(i__k)-  xub(i__k))/  xub(i__k))  <.01) str=TRIM(str)//'k+ '
IF (ABS((xx(i__b)-  xub(i__b))/  xub(i__b))  <.01) str=TRIM(str)//'b+ '
IF (ABS((xx(i__tau)-xub(i__tau))/xub(i__tau))<.01) str=TRIM(str)//'tau+ '
IF (ABS((xx(i__G)-  xub(i__G))/  xub(i__G))  <.01) str=TRIM(str)//'G+ '
IF (ABS((xx(i__pk)- xub(i__pk))/ xub(i__pk)) <.01) str=TRIM(str)//'pk+ '

IF (LEN(TRIM(str))>0) WRITE(unit,'(A60)') 'WARNING: bound hit: '//TRIM(str)

END SUBROUTINE writeboundhit




!****************************************************************************
!
! Purpose:
!  Technical subroutines for csi transformations and input/output procedures
!
! Last changed by Vadym Lepetyuk, October 17, 2004
!
!****************************************************************************
SUBROUTINE readcsi18(filename,csi18)
!DEC$ REAL:8
USE modelsetup
USE chebyshev
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(in) :: filename
TYPE (csi_type), INTENT(inout)::csi18
TYPE (csi_type) ::csif
CHARACTER(LEN=5) :: filesignature
DOUBLE PRECISION, POINTER, DIMENSION(:):: nodesf !DIMENSION(nodes)
INTEGER, PARAMETER :: fileid=32
INTEGER iostatus
INTEGER i

WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'readcsi18:Reading initial csi from a file...'
WRITE(unit1,*) '************************************'

OPEN(unit=fileid,file=filename,status='old',action='read', &
	 form='formatted',position='rewind',iostat=iostatus)
IF (iostatus.NE.0) STOP 'Cannot open csi.txt file'

READ(fileid,*,iostat=iostatus) filesignature
IF (iostatus.NE.0) STOP 'Error reading csi.txt file'
IF (filesignature.NE.'CSI18') STOP 'Wrong csi.txt file type'

READ(fileid,*,iostat=iostatus) csif.nodes
IF (iostatus.NE.0) STOP 'Error reading csi.txt file'
WRITE(unit1,'(A17,I8)') 'readcsi24: nodes', csif.nodes

READ(fileid,*,iostat=iostatus) csif.bmin,csif.bmax,csif.pkmin,csif.pkmax
IF (iostatus.NE.0) STOP 'Error reading csi.txt file'

WRITE(unit1,'(A24,2F10.3)') 'readcsi24: bmin,bmax   ', csif.bmin, csif.bmax
WRITE(unit1,'(A24,2F10.3)') 'readcsi24: pkmin,pkmax ', csif.pkmin, csif.pkmax

CALL initgrid(csif)

ALLOCATE(nodesf(csif.nodes))
READ(fileid,*,iostat=iostatus) nodesf !reading bnodes
IF (iostatus.NE.0) STOP 'Error reading csi.txt file'
IF (MAXVAL(ABS(nodesf-csif.bnodes))>1e-6) STOP 'Nodes in csi.txt are not Chebyshev'

READ(fileid,*,iostat=iostatus) nodesf !reading pknodes
IF (iostatus.NE.0) STOP 'Error reading csi.txt file'
IF (MAXVAL(ABS(nodesf-csif.pknodes))>1e-6) STOP 'Nodes in csi.txt are not Chebyshev'
DEALLOCATE(nodesf)

READ(fileid,*,iostat=iostatus) csif%var
IF (iostatus.NE.0) STOP 'Error reading csi.txt file'
CLOSE(fileid)

DO i=1,n_var
	CALL coef_chebyshev(csif,csif%var(i,:,:),csif%coef(i,:,:))
END DO
CALL convertcsi(csif,csi18)

CALL cleargrid(csif)

END SUBROUTINE readcsi18




!Writing csi into a text file
SUBROUTINE writecsi(csi)
USE modelsetup
IMPLICIT NONE
TYPE (csi_type), INTENT(in)::csi
INTEGER, PARAMETER :: fileid=32
INTEGER iostatus

OPEN(unit=fileid,file='csi.txt',status='replace',action='write', &
	form='formatted',position='rewind',iostat=iostatus)
IF (iostatus.NE.0) STOP 'Cannot write csi.txt file'
WRITE(fileid,'(A5)',iostat=iostatus) 'CSI18'
WRITE(fileid,*,iostat=iostatus) csi%nodes
WRITE(fileid,*,iostat=iostatus) csi%bmin,csi%bmax,csi%pkmin,csi%pkmax
WRITE(fileid,*,iostat=iostatus) csi%bnodes
WRITE(fileid,*,iostat=iostatus) csi%pknodes
WRITE(fileid,*,iostat=iostatus) csi%var
IF (iostatus.NE.0) STOP 'Error writing csi.txt file'
CLOSE(fileid)

END SUBROUTINE writecsi




!Writing csi into a Array Visualizer file
SUBROUTINE writeagl(csi)
USE modelsetup
USE avdef !Defines Array Visualizer's routine interfaces
IMPLICIT NONE
TYPE (csi_type), INTENT(in)::csi
CHARACTER(LEN=3), DIMENSION(i__pk), PARAMETER::&
	name = (/'_cy','__l','gam','__k','__b','tau','__G','_pk'/)
INTEGER iostatus
INTEGER(4) ::aglstatus
INTEGER i
DOUBLE PRECISION, DIMENSION(csi%nodes,csi%nodes)::auxs

DO i=1,i__pk
	auxs=csi%var(i,:,:)
	CALL faglStartWatch (auxs, aglstatus)
	CALL faglSaveAsFile (auxs,'csi_'//name(i)//'.agl', aglstatus)
	CALL faglEndWatch (auxs, aglstatus)
END DO

END SUBROUTINE writeagl





!SUBROUTINE conv2418(csi24,csi18)
!USE modelsetup
!DOUBLE PRECISION, dimension(24,nodes,nodes,nodes), intent(in)::csi24
!DOUBLE PRECISION, dimension(18,nodes,nodes), intent(out)::csi18
!INTEGER, PARAMETER :: & ! Order of the elements in csi24
!    j__cy=1,j__co=2,j__l=3,j__i=4,j__gam=5,&
!    j__k=6,j__b=7,j__tau=8,j__G=9,j__pk=10,&
!    j__cy_G=11,j__cy_pk=12,j__co_G=13,j__co_pk=14,&
!    j__l_G=15,j__l_pk=16,j__i_G=17,j__i_pk=18,&
!    j__k_G=19,j__k_pk=20,j__b_G=21,j__b_pk=22,&
!    j__tau_G=23,j__tau_pk=24
!
!
!csi18(i__cy,:,:)=    csi24(j__cy,4,:,:)
!csi18(i__l,:,:)=     csi24(j__l,4,:,:)
!csi18(i__gam,:,:)=   csi24(j__gam,4,:,:)
!csi18(i__k,:,:)=     csi24(j__k,4,:,:)
!csi18(i__b,:,:)=     csi24(j__b,4,:,:)
!csi18(i__tau,:,:)=   csi24(j__tau,4,:,:)
!csi18(i__G,:,:)=     csi24(j__G,4,:,:)
!csi18(i__pk,:,:)=    csi24(j__pk,4,:,:)
!csi18(i__cy_G,:,:)=  csi24(j__cy_G,4,:,:)
!csi18(i__cy_pk,:,:)= csi24(j__cy_pk,4,:,:)
!csi18(i__l_G,:,:)=   csi24(j__l_G,4,:,:)
!csi18(i__l_pk,:,:)=  csi24(j__l_pk,4,:,:)
!csi18(i__k_G,:,:)=   csi24(j__k_G,4,:,:)
!csi18(i__k_pk,:,:)=  csi24(j__k_pk,4,:,:)
!csi18(i__b_G,:,:)=   csi24(j__b_G,4,:,:)
!csi18(i__b_pk,:,:)=  csi24(j__b_pk,4,:,:)
!csi18(i__tau_G,:,:)= csi24(j__tau_G,4,:,:)
!csi18(i__tau_pk,:,:)=csi24(j__tau_pk,4,:,:)
!
!END SUBROUTINE conv2418



!DOUBLE PRECISION :: gdp
!DOUBLE PRECISION :: cd,A,AFl,AFk,k1n
!!GDP, G/GDP, and gam/GDP at the steady state
!gdp = r*xss(i__k)/(1.+n) + w*xss(i__l)
!WRITE(unit0,'(A10,F10.3)') 'GDP', gdp
!WRITE(unit0,'(A10,F10.3)') 'G/GDP', xss(i__G)/gdp	! Illinois 7%
!WRITE(unit0,'(A10,F10.3)') 'gam/GDP', xss(i__gam)/gdp  ! Illinois 2%

!!Fixed prices vs. Cobb-Douglas prices
!cd=0.3   !The Cobb-Douglas exp. for K in F(K,L)=A*K^cd*L^(1-cd)
!A=10.25  !Procompensctivity
!k1n=xss(i__k)/(1.+n)
!AFl=A*(1.-cd)*(k1n/xss(i__l))**cd
!AFk=A*cd*(xss(i__l)/k1n)**(1.-cd)
!WRITE(unit0,'(A10,2F10.3)') 'Wage', w,AFl
!WRITE(unit0,'(A10,2F10.3)') 'Rent', r,AFk