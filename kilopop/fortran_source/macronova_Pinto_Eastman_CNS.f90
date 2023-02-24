MODULE macronova_Pinto_Eastman_CNS

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
  !   CNS 21.09.2022: remove unsued code    *
  !                                         *
  !   CNS 10.12.2022: update precision in   *
  !         heating rates thermalisation    *
  !         efficiency function.            *
  !                                         *
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

  SUBROUTINE macronova(np, parameters, func_hrate, func_therm, read_hrate, heating_rates_file, Nt, luminosity)

      IMPLICIT NONE
      !--------------------------!
      !-----  I/O parameters ----!
      !--------------------------!
      INTEGER, INTENT(IN) :: np, Nt
      DOUBLE PRECISION, INTENT(IN) :: parameters(np)
      DOUBLE PRECISION, INTENT(OUT) :: luminosity(Nt+1,4)
      LOGICAL, INTENT(IN) :: read_hrate, func_hrate, func_therm
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
      LOGICAL, PARAMETER :: DEBUG = .FALSE.

      !-----------------------!
      !-- file unit numbers --!
      !-----------------------!
      INTEGER, PARAMETER :: uSol=      79
      INTEGER, PARAMETER :: uMod=      80
      INTEGER, PARAMETER :: uLum=      81
      INTEGER, PARAMETER :: uHrate=    83
      INTEGER, PARAMETER :: unit_mac=  85

      INTEGER     ::     it,iostatos,i,j,n,m,nLinesHrate,INFO
      DOUBLE PRECISION :: t0,t1,m_ej,v_med,alpha,e_th,DZ_factor,kappa,Temp0,ye
      DOUBLE PRECISION :: v_max,rho0,E0,F0,tau0,dt,dx,dx2,dx3,L0,hrate
      DOUBLE PRECISION :: qm0,x,y,Lum,tm_p,tm
      DOUBLE PRECISION :: n2,n3,n4,n5,n6,n8,n10,n12,n14,n16,n18,nn1

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
      !--------------------!
      !-- derived values --!
      !--------------------!
      ! [cm/s] maximum velocity, Wollaeger+ (2018), Eq.(12)
      v_max= v_med*128.0d0/63.0d0

      ! [g/cm3] initial central density,Wollaeger+ (2018), Eq.(11)
      rho0=  315.0d0*m_ej/(64.0d0*Pi*(v_max*t0)**3)

      ! [erg/cm3] energy density scale
      E0=    sigma*Temp0**4 * 4.0d0/clight

      ! [erg/s] energy flux scale
      F0= 4.0d0*Pi*clight*v_max*t0*E0/(3.0d0*kappa*rho0)

      ! [s] diffusion timescale
      tau0=  3.0d0*kappa*rho0*(v_max*t0)**2/clight

      ! [t0] timestep size
      dt=    (t1/t0 - 1.0d0)/DBLE(Nt-1)

      ! [R(t)] x-spacing (output only)
      dx=    1.0/DBLE(Nx+1)

      ! auxiliary
      dx2=   dx*dx
      dx3=   dx*dx*dx

      ! [erg/s] lum. scale
      L0= 8.0d0*F0/(1.0d0+2.4d0*dx)

      ! current iteration
      it= 0

      !---------------------!
      !-- allocate arrays --!
      !---------------------!
      ALLOCATE(phi(Nx),phi_p(Nx),rhs(mmax,Nrk),qm(mmax),pm(mmax))
      ALLOCATE(Tnn(Nx),Tnp(0:Nx),Mnn(Nx),Mnp(0:Nx))
      ALLOCATE(AP(Nx*(Nx+1)/2),BP(Nx*(Nx+1)/2),W(Nx),WORK(LWORK))
      ALLOCATE(psi(Nx,Nx),LT(2,Nx),Dm(Nx),Nm2(Nx))

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

        IF(n.GE.1)THEN
           ! Int [x^2(1-x^2)^5 h'n h'n, {x, (n-1)*dx, (n+1)*dx}]
           Tnn(n)= dx*(      (2.0d0/3.0d0+n2*2.0d0) &  ! T_nn
                  + dx2*(     -2.0d0*(1.0d0   +n2*(10.0d0 +n2*5.0d0)) &
                  + dx2*(  20.0d0/7.0d0*(1.0d0 +n2*(21.0d0 +n2*(35.0d0 + n2*7.0d0))) &
                  + dx2*( -20.0d0/9.0d0*(1.0d0   +n2*(36.0d0 +n2*(126.0d0+ n2*(84.0d0+n2*9.0d0))))  &
                  + dx2*(10.0d0/11.0d0+n2*(50.0d0+n2*(300.0d0+n2*(420.0d0+n2*(150.0d0+n2*10.0d0)))) &
                  + dx2*( -2.0d0/13.0d0*(1.0d0   +n2*(78.0d0 +n2*(715.0d0+ n2*(1716.0d0+n2*(1287.0d0&
                  + n2*(286.0d0+n2*13.0d0))))))))))))

           ! Int [x^2(1-x^2)^4 hn hn, {x, (n-1)*dx, (n+1)*dx}]
           ! (Integrate[x^2 (1 - n + x/d)^2 (1 - x^2)^4, {x, (n - 1) d, n d}] +
           !  Integrate[x^2(1-x^2)^4(-x/d+n+1)^2,{x,n d, (n+1)d}] // Simplify)
           !  /. {n->4,d->0.05}
           !  [-]  0.00113393
           Tnn(n)= Tnn(n) + 24.0d0*(&
                   dx3*( 1.0d0/15.0d0*(1.0d0+10.0d0*n2) &
                   + dx2*(-8.0d0/105.0d0*(1.0d0 + 21.0d0*n2 + 35.0d0*n4) &
                   + dx2*((1.0d0/21.0d0 + 12.0d0/7.0d0*n2+6.0d0*n4+4.0d0*n6) &
                   + dx2*(-8.0d0/495.0d0*(1.0d0+55.0d0*n2+330.0d0*n4+462.0d0*n6+165.0d0*n8) &
                   + dx2*( 1.0d0/429.0d0+2.0d0/11.0d0*n2+5.0d0/3.0d0*n4+4.0d0*n6+3.0d0*n8+2.0d0/3.0d0*n10) )))))

           ! Int [x^2(1-x^2)^8 hn hn, {x, (n-1)*dx, (n+1)*dx}]
           ! (Integrate[x^2 (1 - n + x/d)^2 (1 - x^2)^8, {x, (n - 1) d, n d}] +
           !  Integrate[x^2(1-x^2)^8(-x/d+n+1)^2,{x,n d, (n+1)d}] // Simplify)
           ! /. { n->14, d-> 0.05}
           ! [-] 0.000959034
           Mnn(n)= &
                   dx3*(   1.0d0/15.0d0  *(1.0d0+n2*10.0d0) &
                   + dx2*( -16.0d0/105.0d0 *(1.0d0+n2*(21.0d0+n2*35.0d0)) &
                   + dx2*(           2.0d0/9.0d0+n2*(8.0d0+n2*(28.0d0+n2*56.0d0/3.0d0)) &
                   + dx2*(-112.0d0/495.0d0 *(1.0d0+n2*(55.0d0+n2*(330.0d0+n2*(462.0d0+ n2*165.0d0)))) &
                   + dx2*(  70.0d0/429.0d0 *(1.0d0+n2*(78.0d0+n2*(715.0d0+n2*(1716.0d0+ n2*(1287.0d0+286.0d0*n2))))) &
                   + dx2*( -16.0d0/195.0d0 *(1.0d0+n2*(105.0d0+n2*(1365.0d0+n2*(5005.0d0+n2*(6435.0d0+n2*(3003.0d0+n2*455.0d0)))))) &
                   + dx2*(   7.0d0/255.0d0 *(1.0d0+n2*(136.0d0+n2*(2380.0d0+n2*(12376.0d0+n2*(24310.0d0+n2*(19448.0d0+n2*(6188.0d0+n2*680.0d0))))))) &
                   + dx2*(dx2/1995.0d0-16.0d0/2907.0d0       +&
                   n2*(dx2*2.0d0/19.0d0 -16.0d0/2907.0d0*171.0d0  +&
                   n2*(dx2*3.0d0     -16.0d0/2907.0d0*3876.0d0 +&
                   n2*(dx2*136.0d0/5.0d0-16.0d0/2907.0d0*27132.0d0+&
                   n2*(dx2*102.0d0   -16.0d0/2907.0d0*75582.0d0+&
                   n2*(dx2*884.0d0/5.0d0-16.0d0/2907.0d0*92378.0d0+&
                   n2*(dx2*442.0d0/3.0d0-16.0d0/2907.0d0*50388.0d0+&
                   n2*(dx2*408.0d0/7.0d0-16.0d0/2907.0d0*11628.0d0+&
                   n2*(dx2*51.0d0/5.0d0 -16.0d0/2907.0d0*969.0d0 +&
                   n2*(dx2*2.0d0/3.0d0)))))))))&
                   ))))))))

        ENDIF

        ! Tnp(n) is actually T_{n,n+1}
        ! Int [x^2(1-x^2)^5 h'n h'n+1, {x, n*dx, (n+1)*dx}]
        Tnp(n)=   dx*((DBLE(n)**3 - DBLE(1+n)**3)/3.0d0 & ! T_{n,n+1}
                + dx2*((DBLE(1+n)**5 - DBLE(n)**5)    &
                + dx2*((DBLE(n)**7 - DBLE(1+n)**7)*10.0d0/7.0d0 &
                + dx2*((DBLE(1+n)**9 - DBLE(n)**9)*10.0d0/9.0d0 &
                + dx2*((DBLE(n)**11 - DBLE(1+n)**11)*5.0d0/11.0d0 &
                + dx2* (DBLE(1+n)**13 - DBLE(n)**13)/13.0d0)))))

        ! Int [x^2(1-x^2)^4 hn hn+1, {x, n*dx, (n+1)*dx}]
        ! Integrate[x^2*(1-x^2)^4(-x/d+n+1)(x/d-n),{x,n d, (n+1)d}] // Simplify
        ! > /. {n->18, d->0.05}
        ! [-] 3.42611e-6
        Tnp(n)= Tnp(n) + 24.0d0*( &
                dx3*( 1.0d0/60.0d0*(3.0d0+10.0d0*nn1) &
                + dx2*(-2.0d0/105.0d0*(5.0d0+7.0d0*nn1*(4.0d0+5.0d0*nn1)) &
                + dx2*( 1.0d0/84.0d0*(7.0d0+6.0d0*nn1*(9.0d0+7.0d0*nn1*(3.0d0+2.0d0*nn1))) &
                + dx2*(-2.0d0/495.0d0*(9.0d0+11.0d0*nn1*(8.0d0+3.0d0*nn1*(9.0d0+nn1*(12.0d0+5.0d0*nn1)))) &
                + dx2*( 1.0d0/1716.0d0*(11.0d0+13.0d0*nn1*(10.0d0+11.0d0*nn1*(2.0d0+nn1)*(2.0d0+nn1*(3.0d0+2.0d0*nn1)))) &
                ))))))

        ! Int [x^2(1-x^2)^8 hn hn+1, {x, n*dx, (n+1)*dx}]
        ! Integrate[x^2*(1-x^2)^8(-x/d+n+1)(x/d-n),{x,n d, (n+1)d}] // Simplify
        ! > /. {n->18, d->0.05}
        ! [-]  2.13949e-9
        Mnp(n)= dx3*(  1.0d0/60.0d0  *(3.0d0 +10.0d0*nn1) &
                + dx2*( -4.0d0/105.0d0 *(5.0d0 + 7.0d0*nn1*(4.0d0+5.0d0*nn1)) &
                + dx2*(  1.0d0/18.0d0  *(7.0d0 + 6.0d0*nn1*(9.0d0+7.0d0*nn1*(3.0d0+2.0d0*nn1))) &
                + dx2*(-28.0d0/495.0d0 *(9.0d0 +11.0d0*nn1*(8.0d0+3.0d0*nn1*(9.0d0+nn1*(12.0d0+5.0d0*nn1)))) &
                + dx2*( 35.0d0/858.0d0 *(11.0d0+13.0d0*nn1*(10.0d0+11.0d0*nn1*(2.0d0+nn1)*(2.0d0+nn1*(3.0d0+2.0d0*nn1)))) &
                + dx2*( -4.0d0/195.0d0 *(13.0d0+nn1*(180.0d0+13.0d0*nn1*(75.0d0+nn1*(200.0d0+nn1*(270.0d0+7.0d0*nn1*(24.0d0+5.0d0*nn1)))))) &
                + dx2*(  7.0d0/1020.0d0*(15.0d0+34.0d0*nn1*(7.0d0+nn1*(45.0d0+nn1*(150.0d0+nn1*(275.0d0+2.0d0*nn1*(135.0d0+nn1*(63.0d0+10.0d0*nn1))))))) &
                + dx2*( -4.0d0/2907.0d0*(17.0d0+19.0d0*nn1*(16.0d0+17.0d0*nn1*(7.0d0+nn1*(28.0d0+nn1*(65.0d0+nn1*(88.0d0+3.0d0*nn1*(22.0d0+nn1*(8.0d0+nn1)))))))) &
                + dx2*(  1.0d0/7980.0d0*(19.0d0+nn1*(378.0d0+19.0d0*nn1*(168.0d0+nn1*(784.0d0+nn1*(2205.0d0+nn1*&
                 (3822.0d0+nn1*(4004.0d0+nn1*(2376.0d0+7.0d0*nn1*(99.0d0+10.0d0*nn1))))))))))))))))))


      ENDDO !<---------------------------------------------------------- n

      ! Correct corner diagonal elements (imposing boundary conditions)
      Tnn(1)=  Tnn(1) + Tnp(0)
      Mnn(1)=  Mnn(1) + Mnp(0)
      Tnn(Nx)= Tnn(Nx)+ Tnp(Nx)/(1.0d0 + 2.40d0*dx)
      Mnn(Nx)= Mnn(Nx)+ Mnp(Nx)/(1.0d0 + 2.40d0*dx)

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
      AP(1:Nx*(Nx+1)/2)= 0.0d0
      BP(1:Nx*(Nx+1)/2)= 0.0d0
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
      LT(1,1)= 0.0d0
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
         Nm2(m)= 0.0d0
         Dm(m)=  0.0d0
         phi(m)= 0.0d0
         DO n=1,Nx
            x= (DBLE(n)*dx)**2
            n2= DBLE(n*n)
            Nm2(m)= Nm2(m) + n2*(1.0d0-x)**8 * psi(n,m)**2
            Dm(m)=  Dm(m)  + n2*(1.0d0-x)**7 * psi(n,m)
            phi(m)= phi(m) + n2*(1.0d0-x)**8 * psi(n,m)
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
            t_HR(i)= dlog(x*day_in_s)
            HR(i)= dlog(y)
         ENDDO
         CLOSE(uHrate)
      ENDIF

      !-------------------------!
      !-- Main evolution loop --!
      !-------------------------!
      tm=   1.0d0
      tm_p= tm
      phi_p(1:mmax)= phi(1:mmax)

      main_evolution_loop: DO it=0,Nt
         ! --- output t[s] - L[erg/s] - T[K] - R[cm] ------
         Lum= 0.0d0
         DO m=1,mmax
            Lum= Lum + phi(m)*psi(Nx,m)
         ENDDO
         Lum= DABS(Lum * L0)
         x= photospheric_radius(2.0d0/3.0d0,tm,kappa,rho0,t0, v_max)    ! photospheric radius
         y= (Lum/(4.0d0*Pi*sigma))**0.25d0 / DSQRT(x) ! temperature

         IF(Lum > 0.)THEN
            luminosity(it+1,1) = tm*t0
            luminosity(it+1,2) = Lum
            luminosity(it+1,3) = y
            luminosity(it+1,4) = x
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

         IF (func_hrate) THEN
            IF (v_med/clight .GT. 0.5d0) THEN
                v_med = 0.5d0*clight
            END IF
            IF (v_med/clight .LT. 0.05d0) THEN
                v_med = 0.05d0*clight
            END IF
            CALL heating_rate_func(v_med/clight,ye,t0*(tm_p + 0.5d0*dt),hrate)
            IF (func_therm) THEN
                hrate = DZ_factor*calc_t_dep_therm_hrate((t0/day_in_s)*(tm_p + 0.5d0*dt), m_ej, v_max, hrate)
            ELSE
                hrate = e_th*DZ_factor*hrate
            END IF
        ! Look into adding time-dependent thermalisation
         ELSE
            hrate = heating_rate(e_th,alpha,tm_p + 0.5d0*dt, t0, DZ_factor,read_hrate, nLinesHrate,HR, t_HR)
         ENDIF

         DO m=1,mmax
            phi(m)= (phi_p(m) * (1.0d0 - pm(m)*tm_p*0.5d0*dt)  &
                 + dt*qm(m)*(tm_p + 0.5d0*dt) &
                *hrate) &
                 / (1.0d0 + pm(m)*tm*0.5d0*dt)
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

    INTEGER          ivar
    DOUBLE PRECISION eps0,ftot,t_day
    DOUBLE PRECISION:: heating_rate     !< [erg/(g*s)] return value
    DOUBLE PRECISION, PARAMETER :: eps= 1.0d-5
    INTEGER, INTENT(IN):: nLinesHrate
    DOUBLE PRECISION,INTENT(IN):: e_th  ! efficiency
    DOUBLE PRECISION,INTENT(IN):: alpha ! exponent
    DOUBLE PRECISION,INTENT(IN):: t     !< time (in units of t0!)
    DOUBLE PRECISION, INTENT(IN) :: DZ_factor, t_HR(nLinesHrate), t0, HR(nLinesHrate)

    LOGICAL, INTENT(IN) :: read_hrate

    SAVE

    eps0= 2.0d10*(t0/day_in_s)**(-alpha) !< [erg/(g*s)] heating rate scale

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
       IF(DZ_factor > 1.0d0 + eps)THEN
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


  FUNCTION mass_func(x) RESULT(y)

    !***********************************
    !                                  *
    ! Mass function: m(v/vmax>x)/m_ej  *
    !                                  *
    !***********************************

    IMPLICIT NONE
    DOUBLE PRECISION:: y
    DOUBLE PRECISION,INTENT(IN):: x

    y= 16.0d0/315.0d0/(x*x) -x*(1.0d0/3.0d0-3.0d0/5.0d0*x**2+3.0d0/7.0d0*x**4-1.0d0/9.0d0*x**6)

  END FUNCTION mass_func


  FUNCTION calc_t_dep_therm_hrate(time, ejecta_mass, max_ejecta_velocity, hrate) RESULT(eps_tot)

      !************************************************************************
      !                                                                       *
      ! Implement a time-dependent thermalisation based on Oleg's notes       *
      !                                                                       *
      ! 15.12.22 CNS: Changed formatting to declare parameters and dble prec. *
      !                                                                       *
      !************************************************************************

    IMPLICIT NONE


    DOUBLE PRECISION, PARAMETER:: A_alpha = 1.2d-11  ! g cm-3 s
    DOUBLE PRECISION, PARAMETER:: A_beta = 1.3d-11  ! g cm-3 s
    DOUBLE PRECISION, PARAMETER:: A_ff = 0.2d-11  ! g cm-3 s
    ! Assume constant fractions from Wollaeger et al. 2017, which are approximately the time-averaged values from this reference
    DOUBLE PRECISION, PARAMETER:: frac_alpha = 0.05d0
    DOUBLE PRECISION, PARAMETER:: frac_beta = 0.2d0
    DOUBLE PRECISION, PARAMETER:: frac_ff = 0.0d0
    DOUBLE PRECISION, PARAMETER:: frac_gamma = 0.4d0  ! 0.4 or 0.5
    DOUBLE PRECISION, PARAMETER:: kappa_gamma = 0.1d0  ! cm^2/g 0.1 in Oleg from Wollaeger from Barnes 2016

    DOUBLE PRECISION, INTENT(IN):: time, ejecta_mass, max_ejecta_velocity, hrate

    DOUBLE PRECISION:: rho_bar, eta_bar_alpha_sq
    DOUBLE PRECISION:: eta_bar_beta_sq, eta_bar_ff_sq, f_bar_alpha, f_bar_beta, eps_tot
    DOUBLE PRECISION:: tau_bar_gamma, f_bar_gamma, half_vel_t
    DOUBLE PRECISION:: eps_alpha, eps_beta, eps_ff, eps_gamma, f_bar_ff
    DOUBLE PRECISION:: time_s, ejecta_mass_g, max_ejecta_velocity_cms, rho_bar_t
    ! convert inputs to correct units
    time_s = time*day_in_s  ! seconds
    ejecta_mass_g = ejecta_mass*msol  ! grams
    max_ejecta_velocity_cms = max_ejecta_velocity*clight  ! cm/s

    half_vel_t = 0.5d0*max_ejecta_velocity_cms*time_s
    ! g/cm^-3
    rho_bar = 0.14d0*(ejecta_mass_g/((half_vel_t)**3.0d0))
    ! unitless
    tau_bar_gamma =0.035d0*(0.1d0*ejecta_mass_g/((half_vel_t)**2.0d0))

    rho_bar_t = time_s*rho_bar

    ! should be unitless
    eta_bar_alpha_sq = 2.0d0*A_alpha/(rho_bar_t)
    eta_bar_beta_sq = 2.0d0*A_beta/(rho_bar_t)
    eta_bar_ff_sq = 2.0d0*A_ff/(rho_bar_t)

    ! unitless
    f_bar_alpha = dlog(1.0d0 + eta_bar_alpha_sq)/(eta_bar_alpha_sq)
    f_bar_beta = dlog(1.0d0 + eta_bar_beta_sq)/(eta_bar_beta_sq)
    f_bar_ff = dlog(1.0d0 + eta_bar_ff_sq)/(eta_bar_ff_sq)
    f_bar_gamma = 1.0d0 - dexp(-tau_bar_gamma)

    ! determine individual heating rates from each process
    eps_alpha = frac_alpha*hrate
    eps_beta = frac_beta*hrate
    eps_ff = frac_ff*hrate
    eps_gamma = frac_gamma*hrate
    ! combine for total heating rate
    eps_tot = f_bar_alpha*eps_alpha + f_bar_beta*eps_beta + f_bar_ff*eps_ff + f_bar_gamma*eps_gamma

  END FUNCTION calc_t_dep_therm_hrate


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
    IF (y.GT.10.0d0) THEN
       x= 4.0d0/DSQRT(315.0d0*y)
    ELSE
       x= 0.5d0
    ENDIF

    iter: DO i=1,20
       x1= x
       x=  x + (mass_func(x)-y)*x/(2.0d0*mass_func(x)+x*(1.0d0-x*x)**3)
       IF (x.LT.0) x= 4.0d0/DSQRT(315.0d0*y)
       IF (x.GT.1) x= 0.5d0
    ENDDO iter

    x= x*v_max*tm*t0

  END FUNCTION photospheric_radius


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
    DOUBLE PRECISION, PARAMETER :: ftot_fix= 0.5d0


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
       tvec(1)=   0.32464790685222789d0
       tvec(2)=   0.51453225751684539d0
       tvec(3)=   0.81547867224009707d0
       tvec(4)=   1.2924465962305569d0
       tvec(5)=   2.0483898119853476d0
       tvec(6)=   3.2464790685222789d0
       tvec(7)=   5.1453225751684535d0
       tvec(8)=   8.1547867224009707d0
       tvec(9)=   12.924465962305566d0
       tvec(10)=  20.483898119853475d0
       tvec(11)=  32.464790685222788d0
       tvec(12)=  51.453225751684585d0
       tvec(13)=  81.547867224009821d0
       tvec(14)= 129.24465962305572d0

       ! efficiency vector
       fvec(1)=    0.68130602017301889d0
       fvec(2)=    0.67411024354860460d0
       fvec(3)=    0.66603434160559860d0
       fvec(4)=    0.66125114158990717d0
       fvec(5)=    0.64752942779906753d0
       fvec(6)=    0.56476272381257997d0
       fvec(7)=    0.43858356212845145d0
       fvec(8)=    0.33320034385805763d0
       fvec(9)=    0.25248013940862857d0
       fvec(10)=   0.20205909671179934d0
       fvec(11)=   0.16185787817594169d0
       fvec(12)=   0.12420294745382035d0
       fvec(13)=   6.6055034918415953D-002
       fvec(14)=   3.2651573001217751D-002

    !------------------------!
    !-- original data DZ31 --!
    !------------------------!
    CASE(2)
       ! time vector [days]
       tvec(1)=   0.32464790685222789d0
       tvec(2)=   0.51453225751684539d0
       tvec(3)=   0.81547867224009707d0
       tvec(4)=   1.2924465962305569d0
       tvec(5)=   2.0483898119853476d0
       tvec(6)=   3.2464790685222789d0
       tvec(7)=   5.1453225751684535d0
       tvec(8)=   8.1547867224009707d0
       tvec(9)=   12.924465962305566d0
       tvec(10)=  20.483898119853475d0
       tvec(11)=  32.464790685222788d0
       tvec(12)=  51.453225751684585d0
       tvec(13)=  81.547867224009821d0
       tvec(14)= 129.24465962305572d0

       ! efficiency vector
       fvec(1)=    0.72496093096935588d0
       fvec(2)=    0.72970590959526771d0
       fvec(3)=    0.73559449462877358d0
       fvec(4)=    0.74935511733570204d0
       fvec(5)=    0.76287302392115630d0
       fvec(6)=    0.72631878901518698d0
       fvec(7)=    0.66831161898858304d0
       fvec(8)=    0.59647454241251674d0
       fvec(9)=    0.52950421400956926d0
       fvec(10)=   0.46696994547145909d0
       fvec(11)=   0.35804896949618614d0
       fvec(12)=   0.23832932982081834d0
       fvec(13)=   0.13867081443154802d0
       fvec(14)=   4.8199279580718062D-002

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

    IF (yp1 > .99d30) THEN
       y2(1)= 0.0d0
       u(1)=  0.0d0
    ELSE
       y2(1)=-0.5d0
       u(1)=  (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    ENDIF

    DO i=2,n-1
       sig=   (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=     sig*y2(i-1) + 2.0d0
       y2(i)= (sig - 1.0d0)/p
       u(i)=  (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    ENDDO

    IF (ypn > .99d30) THEN
       qn= 0.0d0
       un= 0.0d0
    ELSE
       qn= 0.5d0
       un= (3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    ENDIF

    y2(n)= (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
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
    y= a*ya(klo)+b*ya(khi)+((a**3.0d0-a)*y2a(klo)+(b**3.0d0-b)*y2a(khi))*(h**2)/6.0d0


  END SUBROUTINE splint

END MODULE macronova_Pinto_Eastman_CNS
