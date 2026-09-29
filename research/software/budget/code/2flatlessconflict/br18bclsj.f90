PROGRAM budgetrules
!****************************************************************************
!
!  PROGRAM: budgetrules, version using Least Squares Minimization
!
!  PURPOSE:  main program
!
! Last changed by Vadym Lepetyuk, September 6, 2004
! Copyright by Marco Bassetto, Florin Bidian, Vadym Lepetyuk
! This code can be freely distributed and modified for research purposes only, 
! provided this copyright notice is included in the modified code. 
! Proper credit should be given in all publications arising from
! modifications of this code; this should include a citation of 
! "Politics and Efficiency of Separating Capital and Ordinary Government Budgets"
! by Marco Bassetto with Thomas J. Sargent!
! Note: all routines need to be checked for hidden qk=1, qG=1, qpk=1 assumptions
!
!****************************************************************************
!DEC$ REAL:8
USE modelsetup
USE chebyshev
IMPLICIT NONE
TYPE (csi_type) csi,csix
DOUBLE PRECISION, DIMENSION(n_var) :: xss, xssx
INTEGER iostatus
DOUBLE PRECISION :: u, compens, wedgeG, wedgepk
INTEGER i



!Initialize program logs
unit0=0 !redirect output written to unit0 to screen
unit1=6 !redirect output written to unit1 to file
unit2=6 !redirect output written to unit2 to file
OPEN (unit=unit2,file='logfile.txt',status='replace', &
	action='write',position='rewind',iostat=iostatus)
IF (iostatus.NE.0) WRITE(unit0,*) 'WARNING: Cannot open logfile'

!Initialize Chebyshev interpolation grid
csi%nodes=15  !Number of Chebyshev interpolation nodes
csi%bmin=-3.5; csi%bmax=3.0  !Boundaries for b
csi%pkmin=0.0; csi%pkmax=6.6  !Boundaries for pk
CALL initgrid(csi)

x=1.0; y=0.0 !Initialize x and y
CALL initcsi(csi) !Compute initial guess for csi
!CALL testcsi(csi) !Initial tests (optional)
CALL solvecsi(csi) !Compute csi
CALL solvess(csi,xss) !Compute steady state

IF (unit0.NE.unit1) &
	WRITE(unit0,'(2(A5,F10.3))') 'x=',x, 'y=',y
	CALL writexxs(unit0,'budgetrules: Steady State',n_var,xss)
WRITE(unit0,*) 'Utility', utility(xss,xss)





! Loop over different x
xssx=xss !Keep steady state value (which will be needed for compensation computation)
csix%nodes=csi%nodes
csix%bmin=csi%bmin; csix%bmax=csi%bmax
csix%pkmin=csi%pkmin; csix%pkmax=csi%pkmax
CALL initgrid(csix) !Initialize the additional Chebyshev interpolation grid
csix%var=csi%var !Keep csi values, calculated above
csix%coef=csi%coef !Keep Chebyshev coefficients as well
WRITE(unit0,*) '---------------- budgetrules: Loop over x ----------------'
WRITE(unit0,*) '(x,utility,wedgeG,wedgepk,compensation)'

DO x=1.0,-0.005,-0.005

	CALL adjustboundaries(csi) !Adjust grid if necessary ...
	CALL reinitgrid(csi) !... and recompute Chebyshev interpolation nodes

	CALL convertcsi(csix,csi) !csi from the previous step is our guess for current csi
	!CALL initcsi(csi) !(alternative) instead we may make our own initial guess for csi

	CALL solvecsi(csi)
	CALL solvess(csi,xss)
	u = utility(xss,xss)

	CALL efficiency(xss,wedgeG,wedgepk) 	!Asymptotic Efficiency
	CALL compensation(csi,xssx,compens) 	!Compensation for regime change
	WRITE(unit0,'(5F12.5)') x,u,wedgeG,wedgepk,compens

	csix%bmin=csi%bmin; csix%bmax=csi%bmax
	csix%pkmin=csi%pkmin; csix%pkmax=csi%pkmax
	!csix%bnodes=csi%bnodes; csix%pknodes=csi%pknodes (not necessary)
	!csix%var=csi%var (not necessary)
	csix%coef=csi%coef

END DO

CALL cleargrid(csix)  !Nice but useless step because the program ends anyway
CALL cleargrid(csi)

WRITE(unit1,*) '************************************'
WRITE(unit1,*) 'budgetrules:End of program...'
WRITE(unit1,*) '************************************'

CLOSE(unit2)
CLOSE(unit1)

END PROGRAM budgetrules
