!> @file hratelib.f90
!!
!! module 'hratelib'
!!
!*****************************************************************************
!  Module 'hratelib' for computing heating rates in expanding ejecta         *
!  OK 29.08.2019                                                             *
!  CNS 12.09.2019 modification for inclusion into macronova_Pinto_Eastman_CNS*                                                  *
!*****************************************************************************
MODULE hratelib
IMPLICIT NONE

! grid of velocity and Ye from which approximant is interpolated
DOUBLE PRECISION, PARAMETER :: &
   YE_GRID(10)= (/  .05d0,  .10d0,  .15d0,  .20d0,  .25d0,       &
                    .30d0,  .35d0,  .40d0,  .45d0,  .50d0 /),    &
     V_GRID(6)= (/  .05d0,  .1d0,  .2d0,  .3d0,  .4d0,  .5d0 /)

! approximant coefficients on the grid
DOUBLE PRECISION, PARAMETER, DIMENSION(6, 10) :: &
   E0_GRID= RESHAPE( (/ 5.81d0, 6.76d0, 6.50d0, 6.60d0, 6.60d0, 6.60d0, &
                        6.60d0, 6.60d0, 9.80d0, 9.80d0, 9.80d0, 9.80d0, &
                         6.5d0,    9d0,    1d1,    1d1,    1d1,    1d1, &
                         6.5d0,  6.5d0,  5.8d0, 19.2d0, 18.0d0, 18.0d0, &
                         8.4d0, 29.7d0, 45.2d0, 33.2d0, 33.2d0, 48.0d0, &
                         6.1d0, 18.0d0, 47.1d0, 47.1d0, 74.8d0, 74.8d0, &
                         7.3d0,  7.3d0, 16.3d0, 23.2d0, 43.2d0, 1.50d2, &
                          .1d0, 8.0d-2, 3.0d-2,   .1d0,   .5d0,  .15d0, &
                         1.2d0, 1.4d0,  1.68d0, 1.68d0, 1.68d0, 1.68d0, &
                          .4d0, 1.07d0,  .76d0, 1.42d0, 1.42d0, 1.42d0 /),&
                     (/ 6, 10 /)), &
   ALP_GRID= RESHAPE( (/ 1.35d0, 1.35d0, 1.35d0, 1.35d0, 1.35d0, 1.35d0, &
                         1.35d0, 1.35d0, 1.35d0, 1.35d0, 1.35d0, 1.35d0, &
                         1.32d0, 1.34d0, 1.36d0, 1.32d0, 1.35d0, 1.35d0, &
                         1.36d0, 1.36d0, 1.34d0, 1.34d0, 1.34d0, 1.34d0, &
                         1.368d0,1.34d0, 1.394d0, 1.8d0, 1.66d0, 1.68d0, &
                         1.33d0, 1.33d0, 1.33d0, 1.33d0, 1.374d0,1.374d0,&
                         1.268d0,1.358d0,1.384d0,1.384d0,1.384d0,1.344d0,&
                         1.28d0, 1.28d0, 1.23d0, 1.23d0, 1.1d0,  1.05d0, &
                           8d0,    8d0,  5.8d0,  5.8d0,  5.8d0,  5.8d0,  &
                         1.1d0,  1.1d0,  1.1d0,  1.1d0,  1.1d0,  1.1d0 /), &
                      (/ 6, 10 /)), &
   T0_GRID= RESHAPE( (/1.695d0, 1.34d0, 1.20d0, 1.20d0, 1.20d0,  1.20d0, &
                        1.395d0,1.080d0, .825d0, .825d0, .785d0, .785d0, &
                        1.055d0,  .74d0,  .64d0,  .61d0,  .61d0,  .61d0, &
                          .83d0,  .60d0, .455d0,   .4d0,  .38d0,  .38d0, &
                         .645d0, .355d0,  .22d0,  .17d0,  .14d0,  .10d0, &
                         .540d0,  .29d0,  .18d0,  .13d0, .095d0, .081d0, &
                         .385d0, .235d0,   .1d0,  .06d0, .035d0, .025d0, &
                          .01d0,  .01d0,    0d0, .023d0, .013d0, .013d0, &
                         .195d0, .105d0,  .05d0,  .03d0, .023d0, .021d0, &
                         .128d0, .058d0,  .05d0, .024d0,  .02d0,  .02d0 /), &
                      (/ 6, 10 /)), &
   SIG_GRID=RESHAPE( (/   .1d0,   .1d0,   .1d0,   .1d0,   .1d0,   .1d0, &
                          .1d0,   .1d0,  .08d0,  .08d0, .075d0, .075d0, &
                        .075d0, .075d0,  .09d0,  .06d0,  .08d0,  .08d0, &
                        .075d0, .075d0,  .07d0, .035d0,  .04d0,  .04d0, &
                         .07d0, .021d0,  .02d0, .067d0, .074d0, .074d0, &
                        .096d0, .047d0, .021d0, .021d0, .017d0, .017d0, &
                        .058d0, .094d0, .068d0,  .05d0,  .03d0,  .01d0, &
                        .132d0, .132d0,  .05d0, .004d0,.0015d0,.0015d0, &
                        .144d0, .091d0, .026d0, .019d0, .015d0, .012d0, &
                        .018d0, .004d0, .004d0,.0024d0,.0018d0,.0014d0 /), &
                      (/ 6, 10 /)), &
   ALP1_GRID=RESHAPE((/ 2.98d0, 2.98d0, 2.98d0, 2.98d0, 2.98d0, 2.98d0, &
                        2.98d0, 3.23d0, 3.23d0, 3.23d0, 3.23d0, 3.23d0, &
                        3.23d0, 3.61d0, 3.61d0, 3.97d0, 3.70d0, 3.00d0, &
                        3.23d0, 3.58d0, 3.87d0, 4.00d0, 4.00d0, 4.00d0, &
                        3.23d0,    4d0,    4d0,  3.2d0, 2.41d0, 2.53d0, &
                         3.8d0,  3.8d0,    4d0,    4d0,    4d0,    4d0, &
                         2.4d0,  3.8d0,  3.8d0, 3.21d0, 2.91d0, 3.61d0, &
                         .75d0,  .75d0,  .75d0, 1.46d0, 1.46d0, 1.46d0, &
                           0d0,    0d0,    0d0,    0d0,    0d0,    0d0, &
                         1.4d0,  1.4d0,  1.4d0,    0d0,    0d0,    0d0 /), &
                      (/ 6, 10 /)), &
   T1_GRID= RESHAPE( (/ .205d0, .205d0, .100d0, .080d0, .028d0, .028d0, &
                        .321d0, .103d0, .066d0, .012d0, .012d0, .012d0, &
                        .186d0, .122d0, .052d0, .026d0, .020d0, .020d0, &
                        .186d0, .093d0, .047d0, .037d0, .021d0, .022d0, &
                        .186d0, .131d0, .087d0, .052d0, .043d0, .015d0, &
                         .14d0, .123d0, .089d0, .055d0, .045d0, .031d0, &
                        .264d0,   .1d0,  .07d0, .055d0, .042d0, .033d0, &
                           0d0,    0d0,    0d0, .015d0, .013d0, .013d0, &
                           0d0,    0d0,    0d0,    0d0,    0d0,    0d0, &
                        .051d0, .027d0, .011d0,    0d0,    0d0,    0d0 /), &
                      (/ 6, 10 /)), &
   SIG1_GRID=RESHAPE((/  .16d0,  .18d0,  .10d0,  .20d0, .062d0, .062d0, &
                        .225d0, .061d0, .129d0, .085d0, .085d0, .085d0, &
                        .095d0, .095d0, .034d0, .020d0, .010d0, .010d0, &
                        .095d0, .055d0,  .02d0,  .05d0,  .05d0,  .02d0, &
                        .095d0,  .09d0, .051d0,  .02d0,  .02d0, .015d0, &
                        .055d0, .067d0, .053d0, .032d0, .032d0, .024d0, &
                        .075d0, .044d0,  .03d0,  .02d0,  .02d0, .014d0, &
                          .3d0,   .3d0,   .3d0, .004d0,.0015d0,.0015d0, &
                          .1d0,   .1d0,   .1d0,   .1d0,   .1d0,   .1d0, &
                         .01d0,.0033d0, .002d0,   .1d0,   .1d0,   .1d0 /), &
                      (/ 6, 10 /)), &
   C1_GRID=  RESHAPE((/    0d0,    0d0,    0d0,    0d0,    0d0,    0d0, &
                           0d0,    0d0,    0d0,    0d0,    0d0,    0d0, &
                           0d0,    0d0,    0d0,    0d0,    0d0,    0d0, &
            28.901d0, 28.901d0, 28.138d0, 27.480d0, 27.480d0, 27.480d0, &
            27.909d0, 27.909d0, 27.086d0, 29.888d0, 30.575d0, 29.012d0, &
            27.539d0, 27.539d0, 26.525d0, 26.940d0, 29.289d0,      0d0, &
                 0d0, 26.938d0, 28.036d0, 28.036d0, 27.447d0, 26.663d0, &
            35.576d0, 35.576d0, 35.576d0, 35.459d0, 35.459d0, 34.996d0, &
            20.723d0, 20.723d0, 20.906d0, 20.906d0, 20.906d0, 20.906d0, &
            25.786d0, 28.689d0, 29.012d0, 29.012d0, 29.710d0, 29.710d0 /), &
                      (/ 6, 10 /)), &
   TAU1_GRID=RESHAPE((/    1d0,    1d0,    1d0,    1d0,    1d0,    1d0, &
                           1d0,    1d0,    1d0,    1d0,    1d0,    1d0, &
                           1d0,    1d0,    1d0,    1d0,    1d0,    1d0, &
                        4.15d0, 4.15d0, 5.36d0, 4.75d0, 4.75d0, .432d0, &
                        14.9d0, 11.2d0, 13.2d0, 5.96d0, 3.02d0, 5.35d0, &
                        12.2d0, 12.2d0, 17.2d0, 1.03d0, .613d0,.0086d0, &
                       .0086d0, 13.8d0, 11.4d0, 14.3d0, 13.3d0, 13.3d0, &
                        .038d0, .038d0, .038d0, .038d0, .029d0, .060d0, &
                        1.29d0, 1.29d0, 1.29d0, 1.29d0, 1.29d0, 1.29d0, &
                        12.2d0, 1.38d0, 1.03d0, 1.03d0, 1.03d0, 1.03d0 /), &
                      (/ 6, 10 /)), &
   C2_GRID=  RESHAPE((/    0d0,    0d0,    0d0,    0d0,    0d0,    0d0, &
                           0d0,    0d0,    0d0,    0d0,    0d0,    0d0, &
                           0d0,    0d0,    0d0,    0d0,    0d0,    0d0, &
            23.273d0, 23.445d0, 23.157d0, 22.803d0, 22.803d0, 22.803d0, &
            21.416d0, 22.333d0, 22.481d0, 23.075d0, 24.407d0, 23.673d0, &
                 0d0,      0d0,      0d0, 19.802d0, 22.012d0, 21.045d0, &
                 0d0, 14.138d0, 18.792d0, 19.114d0, 23.810d0, 19.163d0, &
            28.782d0, 28.782d0, 28.324d0, 30.391d0, 30.073d0, 29.428d0, &
            16.760d0, 16.811d0, 16.811d0, 16.811d0, 16.811d0, 16.811d0, &
            23.901d0, 23.901d0, 23.901d0, 23.901d0, 23.901d0, 23.901d0 /), &
                      (/ 6, 10 /)), &
   TAU2_GRID=RESHAPE((/    1d0,    1d0,    1d0,    1d0,    1d0,    1d0, &
                           1d0,    1d0,    1d0,    1d0,    1d0,    1d0, &
                           1d0,    1d0,    1d0,    1d0,    1d0,    1d0, &
                        5.18d0, 4.14d0, 3.71d0, 3.37d0, 3.37d0, 3.37d0, &
                        5.53d0, 6.91d0, 9.15d0, 5.01d0, 4.49d0, 4.49d0, &
                        5.18d0, 5.18d0, 5.18d0, 34.7d0, 8.38d0, 22.6d0, &
                        .864d0, 4.49d0, 95.0d0, 95.0d0, 0.95d0,  146d0, &
                        .043d0, .043d0, .086d0, .025d0, .034d0, .064d0, &
                        3.45d0, 3.45d0, 3.45d0, 3.45d0, 3.45d0, 3.45d0, &
                        25.9d0, 25.9d0, 25.9d0, 25.9d0, 25.9d0, 25.9d0 /), &
                      (/ 6, 10 /))

CONTAINS

  !************************************************************************
  !  Heating rates for arbitrary {v, Ye} at arbitrary time t[s]           *
  !  OK  7.09.2019                                                        *
  !  CNS 15.09.2019                                                       *
  !************************************************************************
  SUBROUTINE heating_rate_func(v, ye, t, h)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN):: v  ! ejecta expansion velocity [c]
  DOUBLE PRECISION, INTENT(IN):: ye ! initial electron fraction
  DOUBLE PRECISION, INTENT(IN):: t  ! time [s]
  DOUBLE PRECISION, INTENT(OUT):: h
  !
  INTEGER:: i1,i2,j1,j2
  DOUBLE PRECISION:: v1, v2, y1, y2, fv, fy, f11, f12, f21, f22
  DOUBLE PRECISION:: e0,alp,t0,sig,alp1,t1,sig1,C1,C2,tau1,tau2
  DOUBLE PRECISION:: a, b, arg1, arg2, val1, val2
  DOUBLE PRECISION, PARAMETER:: oneoverpi = .5d0/acos(0d0)

     IF (v.LT.V_GRID(1).OR.v.GT.V_GRID(SIZE(V_GRID))) STOP "ERROR: v outside the grid"
     find_index_v: DO i1= 1, SIZE(V_GRID)
        IF (v.LT.V_GRID(i1+1).OR.(i1+1).EQ.(SIZE(V_GRID))) EXIT find_index_v
     ENDDO find_index_v
     i2= i1 + 1

     IF (ye.LT.YE_GRID(1).OR.ye.GT.YE_GRID(SIZE(YE_GRID))) STOP "ERROR: Ye outside the grid"
     find_index_ye: DO j1= 1, SIZE(YE_GRID)
        IF (ye.LT.YE_GRID(j1+1).OR.(j1+1).EQ.(SIZE(YE_GRID))) EXIT find_index_ye
     ENDDO find_index_ye
     j2= j1 + 1

     v1= V_GRID(i1)
     v2= V_GRID(i2)
     fv= (v - v1)/(v2 - v1)

     y1= YE_GRID(j1)
     y2= YE_GRID(j2)
     fy= (ye - y1)/(y2 - y1)

     f11= (1d0 - fv)*(1d0 - fy)
     f12= (1d0 - fv)*fy
     f21= fv*(1d0 - fy)
     f22= fv*fy

     e0=   f11*E0_GRID(i1,j1) + f12*E0_GRID(i1,j2) &
         + f21*E0_GRID(i2,j1) + f22*E0_GRID(i2,j2)

     alp=  f11*ALP_GRID(i1,j1) + f12*ALP_GRID(i1,j2) &
         + f21*ALP_GRID(i2,j1) + f22*ALP_GRID(i2,j2)

     t0=   f11*T0_GRID(i1,j1) + f12*T0_GRID(i1,j2) &
         + f21*T0_GRID(i2,j1) + f22*T0_GRID(i2,j2)

     sig=  f11*SIG_GRID(i1,j1) + f12*SIG_GRID(i1,j2) &
         + f21*SIG_GRID(i2,j1) + f22*SIG_GRID(i2,j2)

     alp1= f11*ALP1_GRID(i1,j1) + f12*ALP1_GRID(i1,j2) &
         + f21*ALP1_GRID(i2,j1) + f22*ALP1_GRID(i2,j2)

     t1=   f11*T1_GRID(i1,j1) + f12*T1_GRID(i1,j2) &
         + f21*T1_GRID(i2,j1) + f22*T1_GRID(i2,j2)

     sig1= f11*SIG1_GRID(i1,j1) + f12*SIG1_GRID(i1,j2) &
         + f21*SIG1_GRID(i2,j1) + f22*SIG1_GRID(i2,j2)

     C1=   f11*C1_GRID(i1,j1) + f12*C1_GRID(i1,j2) &
         + f21*C1_GRID(i2,j1) + f22*C1_GRID(i2,j2)

     tau1= f11*TAU1_GRID(i1,j1) + f12*TAU1_GRID(i1,j2) &
         + f21*TAU1_GRID(i2,j1) + f22*TAU1_GRID(i2,j2)

     C2=   f11*C2_GRID(i1,j1) + f12*C2_GRID(i1,j2) &
         + f21*C2_GRID(i2,j1) + f22*C2_GRID(i2,j2)

     tau2= f11*TAU2_GRID(i1,j1) + f12*TAU2_GRID(i1,j2) &
         + f21*TAU2_GRID(i2,j1) + f22*TAU2_GRID(i2,j2)

     a= .5d0 - oneoverpi*datan((t - t0)/sig)
     b= .5d0 + oneoverpi*datan((t - t1)/sig1)

     arg1 = C1 - t/tau1*1d-3
     If (arg1 .LT. -100.0d0) THEN
          val1 = 0.0d0
     ELSE IF (arg1 .GT. 150.0d0) THEN
          val1 = dexp(150.0d0)
     ELSE
          val1 = dexp(arg1)
     END IF

     arg2 = C2 - t/tau2*1d-5
     If (arg2 .LT. -100.0) THEN
          val2 = 0.0d0
     ELSE IF (arg2 .GT. 150.0d0) THEN
          val2 = dexp(150.0d0)
     ELSE
          val2 = dexp(arg2)
     END IF

     h = e0*1e18*(a**alp * b**alp1) + val1 + val2
  END SUBROUTINE heating_rate_func

ENDMODULE hratelib
