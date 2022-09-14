MODULE macronova_Pinto_eastman_CNS

  !******************************************
  !                                         *
  !     MACRONOVA CALCULATOR                *
  !    after Wollaeger+17 (App.A)           *
  !    and Pinto & Eastman'00               *
  !                                         *
  !   original code by Oleg Korobkin        *
  !   12/2017;                              *
  !                                         *
  !   SKR 01/2018: some cleanup +           *
  !                adaptions                *
  !                                         *
  !   SKR 14.01.2018: apply time-dep.       *
  !                   efficiency factors    *
  !                   either FRDM (ivar= 1) *
  !                   or DZ31 (ivar= 2)     *
  !                                         *
  !   CNS 18.06.2018: adapt code for python *
  !                   wrapper to interface  *
  !                   with MNe observations *
  !                   simulation code       *
  !                                         *
  ! ===> read-in heating-rate file MUST BE  *
  !      WITHOUT EFFICIENCY FACTORS         *
  !                                         *
  !******************************************

  !------------------!
  !-- Compilation: --!
  !------------------!
  ! $ gfortran  -O3 -fdefault-real-8 -fdefault-double-8 -ffpe-trap=overflow,underflow,invalid -lm macronova_calculator.f90 -o macronova_calculator.x
  !

  USE physics_constants
  USE hratelib

  CONTAINS

  SUBROUTINE macronova(np, parameters, func_hrate, read_hrate, heating_rates_file, Nt, luminosity)

      IMPLICIT NONE
      !--------------------------!
      !-----  I/O parameters ----!
      !--------------------------!
      INTEGER, INTENT(IN) :: np, Nt
      DOUBLE PRECISION, INTENT(IN) :: parameters(np)
      DOUBLE PRECISION, INTENT(OUT) :: luminosity(Nt+1,4)
      LOGICAL, INTENT(IN) :: read_hrate, func_hrate
      CHARACTER*255, INTENT(IN) :: heating_rates_file

      !--------------------------!
      !-- numerical parameters --!
      !--------------------------!
      INTEGER, PARAMETER :: mmax=  100  ! number of harmonics
      INTEGER, PARAMETER :: Nx=    100  ! x-grid (output only)
      INTEGER, PARAMETER :: LWORK= 3*Nx-1
      INTEGER, PARAMETER :: Nrk=   6    ! number of Runge-Kutta substeps

      !------------------------!
      !-- control parameters --!
      !------------------------!
      ! read precomputed (with quad precision) matrix?
      ! LOGICAL, PARAMETER :: read_matrix = .FALSE.
      LOGICAL, PARAMETER :: DEBUG = .FALSE.

      !-----------------------!
      !-- file unit numbers --!
      !-----------------------!
      INTEGER, PARAMETER :: uSol=      79
      INTEGER, PARAMETER :: uMod=      80
      INTEGER, PARAMETER :: uLum=      81
      INTEGER, PARAMETER :: uHrate=    83
      !INTEGER, PARAMETER :: uMat=      84
      INTEGER, PARAMETER :: unit_mac=  85

      INTEGER     ::     it,iostatos,i,j,n,m,nLinesHrate,INFO
      DOUBLE PRECISION :: t0,t1,m_ej,v_med,alpha,e_th,DZ_factor,kappa,Temp0,ye
      DOUBLE PRECISION :: v_max,rho0,E0,F0,tau0,dt,dx,dx2,dx3,L0,hrate
      DOUBLE PRECISION :: qm0,x,y,Lum,tm_p,tm
      DOUBLE PRECISION :: n2,n3,n4,n5,n6,n8,n10,n12,n14,n16,n18,nn1
      !DOUBLE PRECISION Esol, tau, t_in_days, Flux

      !---------------------------!
      !-- other local variables --!
      !---------------------------!
      ! modes of E(x,t), current and previous timesteps
      DOUBLE PRECISION,ALLOCATABLE::  phi(:),phi_p(:)

      ! auxiliary array for the time integrator
      DOUBLE PRECISION,ALLOCATABLE::  rhs(:,:)

      ! helpers for computing the right-hand side
      DOUBLE PRECISION,ALLOCATABLE::  qm(:),pm(:)

      ! kinetic matrix diagonal and superdiagonal
      DOUBLE PRECISION,ALLOCATABLE::  Tnn(:),Tnp(:)

      ! mass matrix, diagonal / superdiagonal
      DOUBLE PRECISION,ALLOCATABLE::  Mnn(:),Mnp(:)

      ! matrices A and B in Ax=lam*Bx, in packed format
      DOUBLE PRECISION,ALLOCATABLE::  AP(:),BP(:)

      ! eigenvalues and eigenvectors
      DOUBLE PRECISION,ALLOCATABLE::  W(:),psi(:,:)

      ! normalizations and rhs-scalar products with EVs
      DOUBLE PRECISION,ALLOCATABLE::  Nm2(:),Dm(:)

      ! work arrays needed by LAPACK routines
      DOUBLE PRECISION,ALLOCATABLE::  WORK(:),LT(:,:)

      ! heating rates: times and values (read from file)
      DOUBLE PRECISION,ALLOCATABLE::  t_HR(:), HR(:)

      !----------------------------------!
      !-- parse command-line arguments --!
      !----------------------------------!
      ! IF (iargc().NE.1) THEN
      !    PRINT '("Usage: ./macronova_PE.x <PARAMETERS-FILE>")'
      !    PRINT '("example PARAMETERS-FILE: data/MN_model_parameters.dat")'
      !    STOP
      ! ENDIF
      ! CALL getarg(1,parameters_file)

    !******************************************************************************
    !******************************************************************************
    !           Need Modifications Here
    !******************************************************************************
    !******************************************************************************

      !---------------------------------!
      !-- problem-specific parameters --!
      !---------------------------------!
      ! initial time [d]
      t0 = parameters(1)*day_in_s

      ! final time [d]
      t1 = parameters(2)*day_in_s

      ! ejecta mass [Msol]
      m_ej = parameters(3)*msol

      ! median velocity [c]
      v_med = parameters(4)*clight

      ! exponent nuclear heating rate
      alpha = parameters(5)

      ! thermalization factor
      e_th = parameters(6)

      ! heating enhancement-/DZ-factor
      DZ_factor = parameters(7)

      ! effective opacity [cm^2/g]
      kappa = parameters(8)

      ! initial temperature [K]
      Temp0 = parameters(9)

      ! electron fraction, only needed if going to compute the hrate from hratelib
      IF (func_hrate) THEN
          ye = parameters(10)
      ENDIF
    !******************************************************************************
    !******************************************************************************
    !           Need Modifications Here
    !******************************************************************************
    !******************************************************************************


      !-------------------------------------!
      !-- output parameters on the screen --!
      !-------------------------------------!
      ! PRINT '("Simulation parameters:")'
      ! PRINT '(" - initial & final times[d]:        ",2(F9.5,",",F9.5))', &
      !                                              t0/day_in_s,t1/day_in_s
      ! PRINT '(" - ejecta mass [Msun]:              ",F9.5)', m_ej/msol
      ! PRINT '(" - median velocity [c]:             ",F9.5)', v_max/clight
      ! PRINT '(" - gray opacity [cm2/g]:            ",F9.5)', kappa
      ! PRINT '(" - initial temperature [K]:         ",F9.2)', Temp0
      ! PRINT '(" - number of harmonics:             ",I9)', mmax
      ! PRINT '(" - number of timesteps:             ",I9)', Nt
      ! IF (read_hrate) THEN
      !    PRINT '(" - heating rates - from file:        """,A,"""")',&
      !                                            TRIM(heating_rates_file)
      ! ELSE
      !    PRINT '(" - nuclear heating rate:              analytic")'
      !    PRINT '(" - heating rate exponent:           ",F9.5)',alpha
      ! ENDIF
      ! PRINT '(" - thermalization factor:           ",F9.5)', e_th
      ! PRINT '(" - heating enhancement (DZ-factor): ",F9.5)', DZ_factor

      !--------------------!
      !-- derived values --!
      !--------------------!
      ! [cm/s] maximum velocity, Wollaeger+ (2018), Eq.(12)
      v_max= v_med*128./63.

      ! [g/cm3] initial central density,Wollaeger+ (2018), Eq.(11)
      rho0=  315.0*m_ej/(64.0*Pi*(v_max*t0)**3)

      ! [erg/cm3] energy density scale
      E0=    sigma*Temp0**4 * 4/clight

      ! [erg/s] energy flux scale
      F0= 4*Pi*clight*v_max*t0*E0/(3.*kappa*rho0)

      ! [s] diffusion timescale
      tau0=  3.0*kappa*rho0*(v_max*t0)**2/clight

      ! [t0] timestep size
      dt=    (t1/t0 - 1.0)/DBLE(Nt-1)

      ! [R(t)] x-spacing (output only)
      dx=    1.0/DBLE(Nx+1)

      ! auxiliary
      dx2=   dx*dx
      dx3=   dx*dx*dx

      ! [erg/s] lum. scale
      L0= 8.*F0/(1.+2.4*dx)

      ! current iteration
      it= 0

      !---------------------!
      !-- allocate arrays --!
      !---------------------!
      ALLOCATE(phi(Nx),phi_p(Nx),rhs(mmax,Nrk),qm(mmax),pm(mmax))
      ALLOCATE(Tnn(Nx),Tnp(0:Nx),Mnn(Nx),Mnp(0:Nx))
      ALLOCATE(AP(Nx*(Nx+1)/2),BP(Nx*(Nx+1)/2),W(Nx),WORK(LWORK))
      ALLOCATE(psi(Nx,Nx),LT(2,Nx),Dm(Nx),Nm2(Nx))

      !-- open matrix file for reading or writing, depending on the parameter
      ! IF(read_matrix) THEN
      !    Write(*,*) "Debug 1.2"
      !    OPEN (uMat,FILE='matrix.dat',status='OLD')
      !    Write(*,*) "Debug 1.3"
      ! ELSE
      !    Write(*,*) "Debug 1.4"
      !    OPEN (uMat,FILE='matrix.dat',status='UNKNOWN')
      !    Write(*,*) "Debug 1.5"
      ! ENDIF


      DO n=0,Nx
         n2=  DBLE(n*n)
         n3=  DBLE(n)*n2
         n4=  n2*n2
         n5=  DBLE(n)*n4
         n6=  n4*n2
         n8=  n4*n4
         n10= n4*n6
         n12= n6*n6
         n14= n8*n6
         n16= n8*n8
         n18= n8*n10
         nn1= DBLE(n*(n+1))

         ! IF(read_matrix)THEN
         !    READ(uMat,'(4(ES23.16,1X))') Tnn(n),Tnp(n),Mnn(n),Mnp(n)
         ! ELSE
        IF(n.GE.1)THEN
           ! Int [x^2(1-x^2)^5 h'n h'n, {x, (n-1)*dx, (n+1)*dx}]
           Tnn(n)= dx*(      (2./3.+n2*2.) &  ! T_nn
                  + dx2*(     -2.*(1.   +n2*(10. +n2*5.)) &
                  + dx2*(  20./7.*(1.   +n2*(21. +n2*(35. + n2*7.))) &
                  + dx2*( -20./9.*(1.   +n2*(36. +n2*(126.+ n2*(84.+n2*9.))))  &
                  + dx2*(10./11.+n2*(50.+n2*(300.+n2*(420.+n2*(150.+n2*10.)))) &
                  + dx2*( -2./13.*(1.   +n2*(78. +n2*(715.+ n2*(1716.+n2*(1287.&
                  + n2*(286.+n2*13.))))))))))))

           ! Int [x^2(1-x^2)^4 hn hn, {x, (n-1)*dx, (n+1)*dx}]
           ! (Integrate[x^2 (1 - n + x/d)^2 (1 - x^2)^4, {x, (n - 1) d, n d}] +
           !  Integrate[x^2(1-x^2)^4(-x/d+n+1)^2,{x,n d, (n+1)d}] // Simplify)
           !  /. {n->4,d->0.05}
           !  [-]  0.00113393
           Tnn(n)= Tnn(n) + 24.*(&
                   dx3*( 1./15.*(1.+10.*n2) &
                   + dx2*(-8./105.*(1.+21.*n2+35.*n4) &
                   + dx2*((1./21.+12./7.*n2+6.*n4+4.*n6) &
                   + dx2*(-8./495.*(1.+55.*n2+330.*n4+462.*n6+165.*n8) &
                   + dx2*( 1./429.+2./11.*n2+5./3.*n4+4.*n6+3.*n8+2./3.*n10) )))))

           ! Int [x^2(1-x^2)^8 hn hn, {x, (n-1)*dx, (n+1)*dx}]
           ! (Integrate[x^2 (1 - n + x/d)^2 (1 - x^2)^8, {x, (n - 1) d, n d}] +
           !  Integrate[x^2(1-x^2)^8(-x/d+n+1)^2,{x,n d, (n+1)d}] // Simplify)
           ! /. { n->14, d-> 0.05}
           ! [-] 0.000959034
           Mnn(n)= &
                   dx3*(   1./15.  *(1.+n2*10.) &
                   + dx2*( -16./105. *(1.+n2*(21.+n2*35.)) &
                   + dx2*(           2./9.+n2*(8.+n2*(28.+n2*56./3.)) &
                   + dx2*(-112./495. *(1.+n2*(55.+n2*(330+n2*(462.+ n2*165.)))) &
                   + dx2*(  70./429. *(1.+n2*(78.+n2*(715.+n2*(1716.+ n2*(1287.+286.*n2))))) &
                   + dx2*( -16./195. *(1.+n2*(105.+n2*(1365.+n2*(5005.+n2*(6435.+n2*(3003.+n2*455.)))))) &
                   + dx2*(   7./255. *(1.+n2*(136.+n2*(2380.+n2*(12376.+n2*(24310.+n2*(19448.+n2*(6188.+n2*680.))))))) &
                   + dx2*(dx2/1995.-16./2907.       +&
                   n2*(dx2*2./19. -16./2907.*171.  +&
                   n2*(dx2*3.     -16./2907.*3876. +&
                   n2*(dx2*136./5.-16./2907.*27132.+&
                   n2*(dx2*102.   -16./2907.*75582.+&
                   n2*(dx2*884./5.-16./2907.*92378.+&
                   n2*(dx2*442./3.-16./2907.*50388.+&
                   n2*(dx2*408./7.-16./2907.*11628.+&
                   n2*(dx2*51./5. -16./2907.*969.  +&
                   n2*(dx2*2./3.)))))))))&
                   ))))))))

        ENDIF

        ! Tnp(n) is actually T_{n,n+1}
        ! Int [x^2(1-x^2)^5 h'n h'n+1, {x, n*dx, (n+1)*dx}]
        Tnp(n)=   dx*((DBLE(n)**3 - DBLE(1+n)**3)/3. & ! T_{n,n+1}
                + dx2*((DBLE(1+n)**5 - DBLE(n)**5)    &
                + dx2*((DBLE(n)**7 - DBLE(1+n)**7)*10./7. &
                + dx2*((DBLE(1+n)**9 - DBLE(n)**9)*10./9. &
                + dx2*((DBLE(n)**11 - DBLE(1+n)**11)*5./11. &
                + dx2* (DBLE(1+n)**13 - DBLE(n)**13)/13.)))))

        ! Int [x^2(1-x^2)^4 hn hn+1, {x, n*dx, (n+1)*dx}]
        ! Integrate[x^2*(1-x^2)^4(-x/d+n+1)(x/d-n),{x,n d, (n+1)d}] // Simplify
        ! > /. {n->18, d->0.05}
        ! [-] 3.42611e-6
        Tnp(n)= Tnp(n) + 24.0*( &
                dx3*( 1./60*(3+10*nn1) &
                + dx2*(-2./105*(5+7*nn1*(4+5*nn1)) &
                + dx2*( 1./84*(7+6*nn1*(9+7*nn1*(3+2*nn1))) &
                + dx2*(-2./495*(9+11*nn1*(8+3*nn1*(9+nn1*(12+5*nn1)))) &
                + dx2*( 1./1716*(11+13*nn1*(10+11*nn1*(2+nn1)*(2+nn1*(3+2*nn1)))) &
                ))))))

        ! Int [x^2(1-x^2)^8 hn hn+1, {x, n*dx, (n+1)*dx}]
        ! Integrate[x^2*(1-x^2)^8(-x/d+n+1)(x/d-n),{x,n d, (n+1)d}] // Simplify
        ! > /. {n->18, d->0.05}
        ! [-]  2.13949e-9
        Mnp(n)= dx3*(  1./60.  *(3 +10*nn1) &
                + dx2*( -4./105. *(5 + 7*nn1*(4+5*nn1)) &
                + dx2*(  1./18.  *(7 + 6*nn1*(9+7*nn1*(3+2*nn1))) &
                + dx2*(-28./495. *(9 +11*nn1*(8+3*nn1*(9+nn1*(12+5*nn1)))) &
                + dx2*( 35./858. *(11+13*nn1*(10+11*nn1*(2+nn1)*(2+nn1*(3+2*nn1)))) &
                + dx2*( -4./195. *(13+nn1*(180+13*nn1*(75+nn1*(200+nn1*(270+7*nn1*(24+5*nn1)))))) &
                + dx2*(  7./1020.*(15+34*nn1*(7+nn1*(45+nn1*(150+nn1*(275+2*nn1*(135+nn1*(63+10*nn1))))))) &
                + dx2*( -4./2907.*(17+19*nn1*(16+17*nn1*(7+nn1*(28+nn1*(65+nn1*(88+3*nn1*(22+nn1*(8+nn1)))))))) &
                + dx2*(  1./7980.*(19+nn1*(378+19*nn1*(168+nn1*(784+nn1*(2205+nn1*&
                 (3822+nn1*(4004+nn1*(2376+7*nn1*(99+10*nn1))))))))))))))))))

        !WRITE(uMat,'(4(ES23.16,1X))') Tnn(n),Tnp(n),Mnn(n),Mnp(n)
!         ENDIF

    !!!PRINT '(4(ES23.16,1X))', Tnn(n),Tnp(n),Mnn(n),Mnp(n)
      ENDDO !<---------------------------------------------------------- n

      ! Correct corner diagonal elements (imposing boundary conditions)
      Tnn(1)=  Tnn(1) + Tnp(0)
      Mnn(1)=  Mnn(1) + Mnp(0)
      Tnn(Nx)= Tnn(Nx)+ Tnp(Nx)/(1. + 2.4*dx)
      Mnn(Nx)= Mnn(Nx)+ Mnp(Nx)/(1. + 2.4*dx)

      ! Fill in the matrices AP and BP in the packed form:
      ! Two-dimensional storage of the symmetric matrix A:
      !
      !    a11
      !    a21 a22
      !    a31 a32 a33         (aij = aji)
      !    a41 a42 a43 a44
      !
      ! Packed storage of the lower triangle of A:
      ! AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
      !
      ! AP = [ a11, a21, a31, a41, a22, a32, a42, a33, a43, a44 ]
      !
      AP(1:Nx*(Nx+1)/2)= 0.0
      BP(1:Nx*(Nx+1)/2)= 0.0
      DO n=1,Nx
         AP(  n + (n-1)*(2*Nx-n)/2)= Tnn(n)  ! A_{nn}
         BP(  n + (n-1)*(2*Nx-n)/2)= Mnn(n)  ! B_{nn}
         AP(n+1 + (n-1)*(2*Nx-n)/2)= Tnp(n)  ! A_{n+1,n}
         BP(n+1 + (n-1)*(2*Nx-n)/2)= Mnp(n)  ! B_{n+1,n}
      ENDDO

      IF (DEBUG) THEN
         PRINT '("=== Matrix A: ===")'
         DO i=1,Nx
            DO j=1,i
               WRITE (*,'(ES12.5,1X)',ADVANCE='no') AP(i + (j-1)*(2*Nx-j)/2)
            ENDDO
            WRITE (*,*)
         ENDDO
         WRITE (*,*)

         PRINT '("=== Matrix B: ===")'
         DO i=1,Nx
            DO j=1,i
               WRITE (*,'(ES12.5,1X)',ADVANCE='no') BP(i + (j-1)*(2*Nx-j)/2)
            ENDDO
            WRITE (*,*)
         ENDDO
         WRITE (*,*)
      ENDIF

      ! --- solve eigenvalue problem
      ! 1.DPPTRF -> 2.DSPGST -> 3.DSYEV
      !
      CALL DPPTRF('L',Nx,BP,INFO)
      IF (INFO.NE.0) STOP "ERROR: in Cholesky decomposition"

      IF (DEBUG) THEN
         PRINT '("=== Matrix L from LL^T=B: ===")'
         DO i=1,Nx
            DO j=1,i
               WRITE (*,'(ES12.5,1X)',ADVANCE='no') BP(i + (j-1)*(2*Nx-j)/2)
            ENDDO
            WRITE (*,*)
         ENDDO
         WRITE (*,*)
      ENDIF

      ! -2- Prepare eigenvalue problem
      CALL DSPGST (1,'L',Nx,AP,BP,INFO)

      IF (INFO.NE.0) STOP "ERROR: in DSPGST"
      DO i=1,Nx
         DO j=1,i-1
            psi(i,j)= AP(i+(j-1)*(2*Nx-j)/2)
            psi(j,i)= psi(i,j)
         ENDDO
         psi(i,i)= AP(i+(i-1)*(2*Nx-i)/2)
      ENDDO

      IF (DEBUG) THEN
         PRINT '("=== Matrix A~ = (L^-1) A (L^-T) ===")'
         DO i=1,Nx
            DO j=1,Nx
               WRITE (*,'(ES12.5,1X)',ADVANCE='no') psi(i,j)
            ENDDO
            WRITE (*,*)
         ENDDO
         WRITE (*,*)
      ENDIF


      ! -3- Solve eigenvalue problem
      CALL DSYEV ('V','L',Nx,psi,Nx,W,WORK,LWORK,INFO)
      IF (INFO.NE.0) STOP "ERROR: in DSYEV"
      IF (DEBUG) THEN
         PRINT '("=== Eigenvectors psi ===")'
         DO i=1,Nx
            DO j=1,Nx
               WRITE (*,'(ES12.5,1X)',ADVANCE='no') psi(i,j)
            ENDDO
            WRITE (*,*)
         ENDDO
         WRITE (*,*)

         PRINT '("=== Eigenvalues ===")'
         DO i=1,Nx
            WRITE (*,'(ES12.5,1X)',ADVANCE='no') W(i)
         ENDDO
         WRITE (*,*)
      ENDIF

      ! -4- Prepare the LT matrix (B = L*L^T)
      LT(1,1)= 0.0
      LT(2,1)= BP(1)
      DO n=2,Nx
         LT(1,n)= BP(n + (n-2)*(2*Nx-n+1)/2) ! {n,n-1}
         LT(2,n)= BP(n + (n-1)*(2*Nx-n)/2)   ! {n,n}
      ENDDO

      ! -5- Solve LT*z=y
      CALL DTBTRS ('U','N','N',Nx,1,Nx,LT,2,psi,Nx,INFO)
      IF (INFO.NE.0) STOP "ERROR: in DTBTRS"

      IF (DEBUG) THEN
         DO n=2,Nx-1
            x= Tnp(n-1)*psi(n-1,n) + Tnn(n)*psi(n,n) + Tnp(n)*psi(n+1,n)
            y= Mnp(n-1)*psi(n-1,n) + Mnn(n)*psi(n,n) + Mnp(n)*psi(n+1,n)
            print '(101(ES12.5,1X))',W(n),x,W(n)*y
            !print '(101(ES12.5,1X))',W(n),psi(n,1:Nx)
         ENDDO
      ENDIF

      ! --- eigenvalues normalization and projections; initial state
      DO m=1,Nx
         Nm2(m)= 0.0
         Dm(m)=  0.0
         phi(m)= 0.0
         DO n=1,Nx
            x= (DBLE(n)*dx)**2
            n2= DBLE(n*n)
            Nm2(m)= Nm2(m) + n2*(1.-x)**8 * psi(n,m)**2
            Dm(m)=  Dm(m)  + n2*(1.-x)**7 * psi(n,m)
            phi(m)= phi(m) + n2*(1.-x)**8 * psi(n,m)
         ENDDO
         phi(m)= phi(m) / Nm2(m)
         Nm2(m)= Nm2(m) * dx3
         Dm(m)=  Dm(m)  * dx3
         !print '(2(ES12.5,1X))', Nm2(m), Dm(m)
      ENDDO

      ! --- initialize solution and arrays pm and qm in the evolution equation
      qm0 = rho0*t0/E0
      DO m=1,mmax
         qm(m)= qm0 * Dm(m)/Nm2(m)
         pm(m)= W(m)* t0/tau0
      ENDDO

      ! ! --- open files and write headers
      ! PRINT '("Writing output files:")'
      ! PRINT '(" - sol.dat:    snapshots of solution")'
      ! PRINT '(" - modes.dat:  evolution of mode amplitudes")'
      ! PRINT '(" - lum.dat:    time[s], luminosity, temperature, radius")'
      ! PRINT '(" - mac.dat:    time[d], luminosity, temperature, radius, velocity")'
      ! OPEN (uSol,FILE='sol.dat',STATUS='unknown',IOSTAT=iostatos)
      ! IF (iostatos.NE.0) STOP "ERROR: cannot open sol.dat for writing"
      ! WRITE(uSol,'("# Snapshots of solution: 1:x 2:E[erg/cm3] 3:tau 4:total radial flux[erg/s]")')
      ! WRITE(uSol,'("# Simulation parameters:")')
      ! WRITE(uSol,'("# - initial & final times[d]:",2(F9.5,",",F9.5))') t0/day_in_s,t1/day_in_s
      ! WRITE(uSol,'("# - ejecta mass [Msun]:      ",F9.5)') m_ej/msol
      ! WRITE(uSol,'("# - ejecta velocity [c]:     ",F9.5)') v_max/clight
      ! WRITE(uSol,'("# - thermalization factor:   ",F9.5)') e_th
      ! WRITE(uSol,'("# - gray opacity [cm2/g]:    ",F9.5)') kappa
      ! WRITE(uSol,'("# - initial temperature [K]: ",F9.2)') Temp0
      ! WRITE(uSol,'("# - number of harmonics:     ",I9)') mmax
      ! WRITE(uSol,'("# - number of timesteps:     ",I9)') Nt
      !
      ! OPEN (uMod,FILE='modes.dat',STATUS='unknown',IOSTAT=iostatos)
      ! IF (iostatos.NE.0) STOP "ERROR: cannot open modes.dat for writing"
      ! WRITE(uMod,'("# Snapshots of phi_m(t): 1:m 2:phi_m")')
      !
      ! OPEN (uLum,FILE='lum.dat',STATUS='unknown',IOSTAT=iostatos)
      ! IF (iostatos.NE.0) STOP "ERROR: cannot open modes.dat for writing"
      ! WRITE(uLum,'("#  1: time[s]     2: L[erg/s]    3: T [K]       4: R [cm]")')
      !
      ! OPEN(unit_mac,FILE='mac.dat',STATUS='unknown',IOSTAT=iostatos)
      ! IF (iostatos.NE.0) STOP "ERROR: cannot open mac.dat for writing"
      ! WRITE(unit_mac,*)'# 1:time[d]   2:L[erg/s]   3:T[K]   4:R[cm]   5:v[c]'

      ! --- read the heating rates file if requested
      !     Heating rate file format:
      !     - first line: header
      !     - rest of the file: two columns,
      !       in format: 1:time[s] 2:heating rate [erg/g/s]
      !
      ! SKR 14.01.2018: the second column should be the
      !                 NAKED-heating rate, i.e. WITHOUT EFFICIENCY
      !
      IF (read_hrate) THEN
         OPEN(uHrate,FILE=TRIM(heating_rates_file),STATUS='old',IOSTAT=iostatos)
         IF (iostatos.NE.0) THEN
            PRINT '(A,A,A)', "ERROR: cannot open ", &
                  TRIM(heating_rates_file)," for reading"
            STOP
         ENDIF

         ! count number of lines in the file
         nLinesHrate= -1 ! accounts for the header
         read_hrate_file: DO
            READ(uHrate,*,IOSTAT=iostatos)
            IF (iostatos.NE.0) EXIT read_hrate_file
            nLinesHrate= nLinesHrate + 1
         ENDDO read_hrate_file

         ! allocate the arrays and read the file
         ALLOCATE (t_HR(nLinesHrate),HR(nLinesHrate))
         REWIND(uHrate)
         READ(uHrate,*) ! skip header
         DO i=1,nLinesHrate
            READ(uHrate,*) x,y
            t_HR(i)= log(x*day_in_s)
            HR(i)= log(y)
         ENDDO
         CLOSE(uHrate)
      ENDIF

      !-------------------------!
      !-- Main evolution loop --!
      !-------------------------!
      tm=   1.0
      tm_p= tm
      phi_p(1:mmax)= phi(1:mmax)

      main_evolution_loop: DO it=0,Nt

    !      ! --- output E[x,t]
    !      IF (mod(it,(Nt/100)).EQ.0) THEN
    !         WRITE(uSol,'("# t= ",F9.5," [d]")') (tm*t0/day_in_s)
    !         x= 0.0
    !         DO i=0,Nx
    !            Esol= solution_internal_energy(x)
    !            Flux= solution_flux(x)
    !            tau= optical_depth(x,tm)
    !            WRITE(uSol,'(F9.5,3(1X,ES14.7))') x, Esol, tau, Flux
    !
    ! !           Flux=-Flux * E0*clight/(3*kappa*rho0*v_max*t0 * tm**2)
    ! !           tau_optdepth = kappa*rho0*v_max*t0/tm**2*(1-x)
    ! !           WRITE(uSol,'(F9.5,5(1X,ES14.7))') x, Esol, Flux, &
    ! !               x*v_max*t0*tm, & ! current distance from the origin
    ! !               tau_optdepth     ! optical depth
    !            x= x + dx
    !         ENDDO
    !         WRITE(uSol,*);WRITE(uSol,*) ! two empty lines in the end of each block
    !      ENDIF

         ! --- output t[s] - L[erg/s] - T[K] - R[cm] ------
         Lum= 0.0
         DO m=1,mmax
            Lum= Lum + phi(m)*psi(Nx,m)
         ENDDO
         Lum= ABS(Lum * L0)
         x= photospheric_radius(2./3.,tm,kappa,rho0,t0, v_max)    ! photospheric radius
         y= (Lum/(4*Pi*sigma))**0.25 / SQRT(x) ! temperature

         IF(Lum > 0.)THEN
            luminosity(it+1,1) = tm*t0
            luminosity(it+1,2) = Lum
            luminosity(it+1,3) = y
            luminosity(it+1,4) = x
!            WRITE(uLum,'(4(ES14.7,1X))') tm*t0, Lum, y, x
         ENDIF

         ! --- Crank-Nicholson integration step
         !
         !  Differential equation:  y' + pm t y ==  qm t^(1-alpha)
         !
         !  y - yp       (t + dt)*y + t*yp
         !  ------ = -pm -----------------  + qm*(t + dt/2)^(1-alpha)
         !    dt                 2
         !
         tm= tm_p + dt
         !---- Step to add function heating rate call from new hratelib
         IF (func_hrate) THEN
            IF (v_med/clight .GT. 0.5) THEN
                v_med = 0.5*clight
            END IF
            IF (v_med/clight .LT. 0.05) THEN
                v_med = 0.05*clight
            END IF
            CALL heating_rate_func(v_med/clight,ye,t0*(tm_p + 0.5*dt),hrate)
            hrate = e_th * DZ_factor * hrate
        ! Look into adding time-dependent thermalisation
         ELSE
            hrate = heating_rate(e_th,alpha,tm_p + 0.5*dt, t0, DZ_factor,read_hrate, nLinesHrate,HR, t_HR)
         ENDIF

         DO m=1,mmax
            phi(m)= (phi_p(m) * (1.0 - pm(m)*tm_p*0.5*dt)  &
                 + dt*qm(m)*(tm_p + 0.5*dt) &
                *hrate) &
                 / (1.0 + pm(m)*tm*0.5*dt)
         ENDDO
         tm_p= tm
         phi_p(1:mmax)= phi(1:mmax)

      ENDDO main_evolution_loop

      !-------------!
      !-- cleanup --!
      !-------------!
      !---------------------!
      !-- deallocate arrays --!
      !---------------------!
      DEALLOCATE(phi,phi_p,rhs,qm,pm)
      DEALLOCATE(Tnn,Tnp,Mnn,Mnp)
      DEALLOCATE(AP,BP,W,WORK)
      DEALLOCATE(psi,LT,Dm,Nm2)
      IF (read_hrate) THEN
        DEALLOCATE(t_HR,HR)
      ENDIF



  END SUBROUTINE macronova

  !
  ! This function returns either analytic or interpolated heating rate
  ! Uses parameters:
  ! - read_hrate
  ! - e_th        : thermalization constant
  ! - t0          : [s] reference time
  ! - day_in_s    : [s] one day
  ! - alpha       : rate decay exponent
  ! - t_HR, HR    : arrays with discrete times and heating rates
  ! - nLinesHrate : number of entries in the heating rates file
  ! - DZ_factor   : mimicking DZ31 heating
  !

  FUNCTION heating_rate(e_th,alpha,t, t0, DZ_factor,read_hrate, nLinesHrate,HR, t_HR)

    !*************************
    !                        *
    ! NET heating rate       *
    !                        *
    !*************************

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: eps= 1.D-5
    INTEGER, INTENT(IN):: nLinesHrate
    DOUBLE PRECISION:: heating_rate     !< [erg/(g*s)] return value
    DOUBLE PRECISION,INTENT(IN):: e_th  ! efficiency
    DOUBLE PRECISION,INTENT(IN):: alpha ! exponent
    DOUBLE PRECISION,INTENT(IN):: t     !< time (in units of t0!)
    DOUBLE PRECISION, INTENT(IN) :: DZ_factor, t_HR(nLinesHrate), t0, HR(nLinesHrate)

    LOGICAL, INTENT(IN) :: read_hrate

    SAVE

    INTEGER          ivar
    DOUBLE PRECISION eps0,ftot,t_day

    eps0= 2.e10*(t0/day_in_s)**(-alpha) !< [erg/(g*s)] heating rate scale

    IF (read_hrate) THEN
       heating_rate= exp(linear_interp(log(t*t0),t_HR,HR,nLinesHrate))
    ELSE
       heating_rate= eps0*t**(-alpha)
    ENDIF
    heating_rate= heating_rate*DZ_factor

    !-----------------------------!
    !-- apply efficiency factor --!
    !-----------------------------!
    IF(.NOT.read_hrate)THEN
       ftot= e_th ! the one from input file
    ELSE
       IF(DZ_factor > 1.0 + eps)THEN
          ivar= 2
       ELSE
          ivar= 1
       ENDIF
       t_day= t*t0/day_in_s
       CALL av_efficiency_FRDM_or_DZ(ivar,t_day,ftot)
    ENDIF

    heating_rate= ftot*heating_rate

  END FUNCTION heating_rate


  FUNCTION linear_interp(x,x_arr,y_arr,N) RESULT(y)

    !*********************************
    !                                *
    ! Linear interpolation function  *
    !                                *
    !*********************************

    IMPLICIT NONE

    DOUBLE PRECISION:: y
    INTEGER,INTENT(IN):: N
    DOUBLE PRECISION,INTENT(IN):: x,x_arr(N),y_arr(N)
    INTEGER:: k

    k= get_index(x,x_arr,N)
    y= y_arr(k) + (y_arr(k+1) - y_arr(k))*(x - x_arr(k)) &
         / (x_arr(k+1) - x_arr(k))

  END FUNCTION linear_interp


  FUNCTION get_index(x,A,N) RESULT(k)

    !*********************************************
    !                                            *
    ! Find index i such that A(i) <= x < A(i+1)  *
    ! A(1:N) is a non-decreasing array           *
    !                                            *
    !*********************************************

    IMPLICIT NONE

    INTEGER:: k,N
    DOUBLE PRECISION,INTENT(IN):: x, A(N)

    DO k=N-1,2,-1
       IF (x.GE.A(k)) EXIT
    ENDDO

  END FUNCTION get_index


  FUNCTION optical_depth(x,tm, kappa, rho0, v_max, t0) RESULT(y)

    !************************************
    !                                   *
    ! Optical depth as a function of x  *
    !                                   *
    !************************************

    IMPLICIT NONE
    DOUBLE PRECISION:: y
    DOUBLE PRECISION,INTENT(IN):: x
    DOUBLE PRECISION,INTENT(IN):: tm
    DOUBLE PRECISION, INTENT(IN) :: kappa, rho0, v_max, t0

    y= x*(x*x-1.) - 3./5.*(x**5-1.) + (x**7 - 1.)/7.
    y= y * kappa*rho0*v_max*t0/(tm*tm)

  END FUNCTION optical_depth


  FUNCTION mass_func(x) RESULT(y)

    !***********************************
    !                                  *
    ! Mass function: m(v/vmax>x)/m_ej  *
    !                                  *
    !***********************************

    IMPLICIT NONE
    DOUBLE PRECISION:: y
    DOUBLE PRECISION,INTENT(IN):: x

    y= 16./315./(x*x) -x*(1./3.-3./5.*x**2+3./7.*x**4-1./9.*x**6)

  END FUNCTION mass_func


  FUNCTION photospheric_radius(tau,tm,kappa,rho0,t0, v_max) RESULT(x)

    !****************************
    !                           *
    ! Recover normalized radius *
    !                           *
    !****************************

    IMPLICIT NONE
    DOUBLE PRECISION:: x,x1,y
    DOUBLE PRECISION,INTENT(IN):: tau !< photospheric optical depth
    DOUBLE PRECISION,INTENT(IN):: tm  !< normalized time: t/t0
    DOUBLE PRECISION, INTENT(IN):: kappa, rho0, t0, v_max
    INTEGER:: i

    y= tau * (tm*tm) / (rho0*kappa*v_max*t0)
    IF (y.GT.10.0) THEN
       x= 4./SQRT(315.*y)
    ELSE
       x= 0.5
    ENDIF

    iter: DO i=1,20
       x1= x
       x=  x + (mass_func(x)-y)*x/(2.*mass_func(x)+x*(1.-x*x)**3)
       IF (x.LT.0) x= 4./SQRT(315.*y)
       IF (x.GT.1) x= 0.5
    ENDDO iter

    x= x*v_max*tm*t0

  END FUNCTION photospheric_radius


  FUNCTION solution_internal_energy (x,Nx,N,phi, psi, E0) RESULT(y)

    !*********************************************
    !                                            *
    ! Assemble solution at radial coordinate x_i *
    !                                            *
    !*********************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nx, N
    DOUBLE PRECISION y
    DOUBLE PRECISION, INTENT(IN):: x        !< grid coordinate
    DOUBLE PRECISION, INTENT(IN):: phi(Nx), psi(Nx,Nx), E0

    INTEGER i
    i= MAX(INT((x+1e-12)*DBLE(Nx+1)),1)
    y= E0 * SUM( psi(i ,1:N) * phi(1:Nx)) * (1.-x*x)**4

  END FUNCTION solution_internal_energy


  FUNCTION solution_flux (x,dx,F0,Nx,phi,psi) RESULT(y)

    !*********************************************
    !                                            *
    ! Assemble solution at radial coordinate x_i *
    !                                            *
    !*********************************************

    IMPLICIT NONE

    DOUBLE PRECISION y
    DOUBLE PRECISION, INTENT(IN):: x,dx,F0  !< grid coordinate
    INTEGER, INTENT(IN) :: Nx
    DOUBLE PRECISION, INTENT(IN):: psi(Nx,Nx),phi(Nx)

    INTEGER i,m
    DOUBLE PRECISION z, dzdx
    i= MAX(INT((x+1e-12)*DBLE(Nx+1)),1)
    y= 0.0
    DO m=1,Nx
       z= psi(i,m)
       IF(m.EQ.1) THEN
          dzdx= (psi(2,m) - psi(1,m))/dx
       ELSEIF(m.EQ.Nx) THEN
          dzdx= (psi(Nx,m)/(1.+2.4*dx) - psi(Nx-1,m))/(2.*dx)
       ELSE
          dzdx= (psi(i+1,m) - psi(i-1,m))/(2.*dx)
       ENDIF
       y= y + phi(m)*(8.*x*z - (1.-x*x)*dzdx)
    ENDDO
    y= y * F0 * x*x

  END FUNCTION solution_flux


  SUBROUTINE av_efficiency_FRDM_or_DZ(ivar,t_d,f_tot)

    !*********************************************************
    !                                                        *
    ! use a spline fit to the average heating efficiency;    *
    ! data for fit are an average over the five NSNS cases   *
    ! shown in Rosswog et al. (2017);                        *
    !                                                        *
    ! SKR 14.01.2018:                                        *
    !                 a) separate efficiencies FRDM and DZ31 *
    !                 b) do NOT apply DZ-factor here         *
    !                                                        *
    ! ivar decides:   ivar= 0 ==> constant efficiency        *
    !                 ivar= 1 ==> FRDM-efficieny, t-dep.     *
    !                 ivar= 2 ==> DZ31-efficieny, t-dep.     *
    !                                                        *
    !*********************************************************

    IMPLICIT NONE

    INTEGER,          INTENT(IN) :: ivar  ! const/FRDM/DZ31 ?
    DOUBLE PRECISION, INTENT(IN) :: t_d   ! time [days]
    DOUBLE PRECISION, INTENT(OUT):: f_tot ! fitted average efficiency

    INTEGER,          PARAMETER :: nvec= 14
    DOUBLE PRECISION, PARAMETER :: ftot_fix= 0.5


    DOUBLE PRECISION tvec(nvec),fvec(nvec),y2(nvec)
    DOUBLE PRECISION yp1,ypn

    !----------------------------!
    !-- time dependent or not? --!
    !----------------------------!
    SELECT CASE(ivar)
    CASE(0)
       f_tot= ftot_fix

    !-----------------------!
    !-- original data FRDM--!
    !-----------------------!
    CASE(1)
       ! time vector [days]
       tvec(1)=   0.32464790685222789
       tvec(2)=   0.51453225751684539
       tvec(3)=   0.81547867224009707
       tvec(4)=   1.2924465962305569
       tvec(5)=   2.0483898119853476
       tvec(6)=   3.2464790685222789
       tvec(7)=   5.1453225751684535
       tvec(8)=   8.1547867224009707
       tvec(9)=   12.924465962305566
       tvec(10)=  20.483898119853475
       tvec(11)=  32.464790685222788
       tvec(12)=  51.453225751684585
       tvec(13)=  81.547867224009821
       tvec(14)= 129.24465962305572

       ! efficiency vector
       fvec(1)=    0.68130602017301889
       fvec(2)=    0.67411024354860460
       fvec(3)=    0.66603434160559860
       fvec(4)=    0.66125114158990717
       fvec(5)=    0.64752942779906753
       fvec(6)=    0.56476272381257997
       fvec(7)=    0.43858356212845145
       fvec(8)=    0.33320034385805763
       fvec(9)=    0.25248013940862857
       fvec(10)=   0.20205909671179934
       fvec(11)=   0.16185787817594169
       fvec(12)=   0.12420294745382035
       fvec(13)=   6.6055034918415953E-002
       fvec(14)=   3.2651573001217751E-002

    !------------------------!
    !-- original data DZ31 --!
    !------------------------!
    CASE(2)
       ! time vector [days]
       tvec(1)=   0.32464790685222789
       tvec(2)=   0.51453225751684539
       tvec(3)=   0.81547867224009707
       tvec(4)=   1.2924465962305569
       tvec(5)=   2.0483898119853476
       tvec(6)=   3.2464790685222789
       tvec(7)=   5.1453225751684535
       tvec(8)=   8.1547867224009707
       tvec(9)=   12.924465962305566
       tvec(10)=  20.483898119853475
       tvec(11)=  32.464790685222788
       tvec(12)=  51.453225751684585
       tvec(13)=  81.547867224009821
       tvec(14)= 129.24465962305572

       ! efficiency vector
       fvec(1)=    0.72496093096935588
       fvec(2)=    0.72970590959526771
       fvec(3)=    0.73559449462877358
       fvec(4)=    0.74935511733570204
       fvec(5)=    0.76287302392115630
       fvec(6)=    0.72631878901518698
       fvec(7)=    0.66831161898858304
       fvec(8)=    0.59647454241251674
       fvec(9)=    0.52950421400956926
       fvec(10)=   0.46696994547145909
       fvec(11)=   0.35804896949618614
       fvec(12)=   0.23832932982081834
       fvec(13)=   0.13867081443154802
       fvec(14)=   4.8199279580718062E-002

    END SELECT

    ! for t-dep. efficiencies...
    IF(ivar > 0)THEN

       !-----------------------------------!
       !-- calculate spline coefficients --!
       !-----------------------------------!
       ! only needs to be calculated once, but do not bother here...
       yp1= 1.D30
       ypn= 1.D30
       CALL spline(tvec,fvec,nvec,yp1,ypn,y2)

       !-----------------!
       !-- interpolate --!
       !-----------------!
       CALL splint(tvec,fvec,y2,nvec,t_d,f_tot)

       ! keep constant instead of extrapolating
       IF(t_d > tvec(nvec))f_tot= fvec(nvec)

    ENDIF

  END SUBROUTINE av_efficiency_FRDM_or_DZ


  SUBROUTINE spline(x,y,n,yp1,ypn,y2)

    !*******************************************************************
    !                                                                  *
    ! Subroutine adapted from "Numerical Recices"; SKR 02.05.2015      *
    !                                                                  *
    ! Given arrays x(1:n) and y(1:n) containing a tabulated function,  *
    ! i.e., yi = f(xi), with x1 < x2 < ... < xN, and given values yp1  *
    ! and ypn for the first derivative of the interpolating function at*
    ! points 1 and n, respectively, this routine returns array y2(1:n) *
    ! of length n which contains the second derivatives of the         *
    ! interpolating function at the tabulated points xi. If yp1 and/or *
    ! ypn are equal to 10^30 or larger, the routine is signaled to set *
    ! the corresponding boundary condition for a natural spline, with  *
    ! zero second derivative on that boundary. Parameter: NMAX is the  *
    ! largest anticipated value of n.                                  *
    !                                                                  *
    ! Please notice that the program spline is called only once to     *
    ! process an entire tabulated function in arrays xi and yi. Once   *
    ! done, the values of the interpolated function for any value of x *
    ! are obtained by calls (as many as desired) to a separate         *
    ! routine splint.                                                  *
    !                                                                  *
    !*******************************************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: n
    DOUBLE PRECISION, INTENT(IN)  :: x(n),y(n),yp1,ypn
    DOUBLE PRECISION, INTENT(OUT) :: y2(n)

    INTEGER, PARAMETER            :: NMAX=500

    INTEGER                       :: i,k
    DOUBLE PRECISION              :: p,qn,sig,un,u(NMAX)

    IF (yp1 > .99e30) THEN
       y2(1)= 0.
       u(1)=  0.
    ELSE
       y2(1)=-0.5
       u(1)=  (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    ENDIF

    DO i=2,n-1
       sig=   (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=     sig*y2(i-1) + 2.
       y2(i)= (sig - 1.)/p
       u(i)=  (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    ENDDO

    IF (ypn > .99e30) THEN
       qn= 0.
       un= 0.
    ELSE
       qn= 0.5
       un= (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    ENDIF

    y2(n)= (un-qn*u(n-1))/(qn*y2(n-1)+1.)
    DO k=n-1,1,-1
       y2(k)= y2(k)*y2(k+1) + u(k)
    ENDDO

  END SUBROUTINE spline


  SUBROUTINE splint(xa,ya,y2a,n,x,y)

    !***************************************
    !                                      *
    ! Subroutine from "Numerical Recices", *
    ! adapted as by SKR 16.01.2013         *
    !                                      *
    ! Given the arrays xa(1:n) and ya(1:n) *
    ! of length n, which tabulate a        *
    ! function (with the xai 's in order), *
    ! and given the array y2a(1:n), which  *
    ! is the output from  the "spline"     *
    ! subroutine, and given a value of x,  *
    ! this routine returns a cubic-spline  *
    ! interpolated value y at x.           *
    !                                      *
    !***************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN)           :: n
    DOUBLE PRECISION, INTENT(IN)  :: xa(n),ya(n),y2a(n),x
    DOUBLE PRECISION, INTENT(OUT) :: y

    INTEGER                       :: k,khi,klo
    DOUBLE PRECISION              :: a,b,h

    klo= 1
    khi= n

    ! figure out bracketing indices
1   IF (khi-klo > 1) THEN
       k= (khi + klo)/2
       IF(xa(k) > x)THEN
          khi= k
       ELSE
          klo= k
       ENDIF
       goto 1
    ENDIF

    h= xa(khi)-xa(klo)

    IF (h == 0.) STOP 'bad xa input in splint'

    a= (xa(khi)-x)/h
    b= (x-xa(klo))/h
    y= a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.


  END SUBROUTINE splint

END MODULE macronova_Pinto_eastman_CNS
