!------------------------------------------------------------------------------
! EPFL/Swiss Plasma Center
!------------------------------------------------------------------------------
!
! MODULE: IIEE 
!
!> @author
!> S. Guinchard - EPFL/SPC
!
! DESCRIPTION:
!> Module handling ion induced electron emissions (IIEE)
!> following Schou's model (for Kinetic emissions) and Auger neutralisation 
!> for potential emissions 
!------------------------------------------------------------------------------

MODULE iiee

        USE particletypes
        USE constants
        USE basic
        USE random_distr 
        !USE factorial      ASK GUILLAUME HOW TO COMBINE MODULES 
        !#include "mkl_vsl.f90"

        IMPLICIT NONE
     !> All the coefficients below were obtained by piecewise
     !> fit of dE/dx curve for stainless steel
        REAL(KIND = db), DIMENSION(2) :: coefficients_1H  = (/0.0, 2.9778E02 /)
        REAL(KIND = db), DIMENSION(2) :: coefficients_1He = (/0.1273010048, 1.70478995200000E02 /)
        REAL(KIND = db), DIMENSION(2) :: coefficients_1Ne = (/0.0518524160, 2.45927583999999E02 /)
        REAL(KIND = db), DIMENSION(4) :: coefficients_2 = (/1.834520315818557E02,1.320304216355084E05, &
                                                               -8.583700602370013E06,3.526140145560557E08/)
        REAL(KIND = db), DIMENSION(4) :: coefficients_3 = (/2.471679999999999E02,9.695466666666670E04, &
                                                                -2.475200000000003E06, 2.325333333333340E07/)
        REAL(KIND = db), DIMENSION(4) :: coefficients_4 = (/2.533904454349683E03,-1.766016382825937E05,&
                                                                8.202640024592019E06, -1.125320217235288E08/)
        REAL(KIND = db), DIMENSION(4) :: coefficients_5 = (/8.786057142856745E02,2.856595238095524e+04, &
                                                              -1.834285714286391E05, 3.333333333338620E05/)
        REAL(KIND = db) :: yield

        ! INTEGER indpion, indpelec !< indices of ions and electrons respectively
        ! the above integers would not be used providing that an index iiee_id
        ! is defined

CONTAINS
        
   !---------------------------------------------------------------------------
   !> @author
   !> Salomon Guinchard EPFL/SPC
   !
   ! DESCRIPTION
   !> function to determine the number of electrons
   !> to add to a given species as a function of th number
   !> of lost ions
   !
   !--------------------------------------------------------------------------
SUBROUTINE ion_induced(pion, losthole, pelec, nblostparts)
    
    TYPE(particles), INTENT(INOUT):: pion, pelec !< ion and electrons parts
    REAL(KIND = db), DIMENSION(3) :: last_pos    !< last position for lost ion (revert push)
    INTEGER, DIMENSION(pion%Nploc):: losthole    !< indices of lost ions
    INTEGER :: j, nblostparts, Nploc             !< loop index and #lost particles
    INTEGER :: parts_size_increase, nbadded  
    REAL(KIND = db) :: dr                        !< dr displacement for created parts
                                                 !< will be changed to normal disp                       
   ! nblostparts=2 !< #lost particles (set to 2 for test purposes)
    nbadded=nblostparts !< # particles to add to electron species
    dr = 1E-3/rnorm
    
    IF(pelec%Nploc + nbadded .gt. size(pelec%Z,1)) THEN
      parts_size_increase=Max(floor(0.1*size(pelec%Z,1)),nbadded)
      CALL change_parts_allocation(pelec, parts_size_increase)
    END IF

    
    DO j=1,nblostparts 
        pelec%Nploc = pelec%Nploc+1
        Nploc = pelec%Nploc
        last_pos = revert_push(pion, losthole(j))

      !  pelec%R(Nploc)    = pion%R(losthole(j)) + dr
      !  pelec%Z(Nploc)    = pion%Z(losthole(j)) 
      !  pelec%THET(Nploc) = pion%THET(losthole(j))

        pelec%R(Nploc)     = last_pos(1)
        pelec%THET(Nploc)  = last_pos(2)
        pelec%Z(Nploc)     = last_pos(3)
        pelec%UR(Nploc)    = 0
        pelec%UZ(Nploc)    = 0
        pelec%UTHET(NPLOC) = 0


    END DO 
    

END SUBROUTINE ion_induced


   !---------------------------------------------------------------------------
   !> @author
   !> Salomon Guinchard EPFL/SPC
   !
   ! DESCRIPTION
   !> reverts Buneman algorithm over one time step
   !> to obtain ion position right before recapture
   !> and hence new electron positions
   !
   !--------------------------------------------------------------------------
FUNCTION revert_push(pion, partid)
    USE basic, ONLY: dt, tnorm 
    REAL(KIND=db), DIMENSION(3)::  revert_push   
    TYPE(particles), INTENT(INOUT):: pion !> species: ions
    INTEGER :: partid !> id of particle to reverse position               
    ! We should try to reverse the angle
    ! else, one simple and hence maybe temporary
    ! method is to reverse to previous pos using UR/UTHET*dt
    
    revert_push(1)  = pion%R(partid) - pion%UR(partid)*dt 
    revert_push(2)  = pion%THET(partid) -1/pion%R(partid)* pion%UTHET(partid)*dt
    revert_push(3)  = pion%Z(partid) -pion%UZ(partid)*dt 

    ! BELOW WE TRY TO REVERSE THE ANGLE 
    ! REMAINS TO BE DONE 


END FUNCTION revert_push 


   !---------------------------------------------------------------------------
   !> @author
   !> Salomon Guinchard EPFL/SPC
   !
   ! DESCRIPTION
   !> Evaluate a polynomial at a given point 
   !> with its coefficients provided in an array 
   !> s.t lowest order coeff = 1st element 
   !
   !--------------------------------------------------------------------------
REAL(KIND = db) FUNCTION eval_polynomial(coefficients, value)
    REAL(KIND = db), DIMENSION(:) :: coefficients !< polynomial (e.g fitted yield) coeffs
    REAL(KIND = db) :: value                      !< point where to evaluate polyn
    INTEGER :: ii

    eval_polynomial = 0
    DO ii=1, size(coefficients)
      eval_polynomial = eval_polynomial+coefficients(ii)*value**(ii-1)
    END DO 
END FUNCTION eval_polynomial


   !---------------------------------------------------------------------------
   !> @author
   !> Salomon Guinchard EPFL/SPC
   !
   ! DESCRIPTION
   !> Gives the theoretical value for the electron yield 
   !> as a function of the energy of the incident ion and
   !> the type of neutral gas  
   !
   !--------------------------------------------------------------------------
REAL(KIND = db) FUNCTION compute_yield(energy, neuttype_id)
    REAL(KIND = db) :: energy
    INTEGER :: neuttype_id
    
    IF(energy.le. 1E-3 ) THEN 
            SELECT CASE(neuttype_id)
              CASE(1)
                compute_yield = eval_polynomial(coefficients_1H, energy)
              CASE(2)
                compute_yield = eval_polynomial(coefficients_1He, energy)
              CASE(3)
                compute_yield = eval_polynomial(coefficients_1Ne, energy)
              CASE DEFAULT
                compute_yield = eval_polynomial(coefficients_1H, energy)
       
            END SELECT
        ELSE IF(energy.gt. 1E-3 .and. energy.le. 1E-2) THEN
                compute_yield = eval_polynomial(coefficients_2,energy)

        ELSE IF(energy.gt. 1E-2  .and. energy.le. 2E-2 ) THEN
                compute_yield = eval_polynomial(coefficients_3,energy)

        ELSE IF(energy.gt. 2E-2 .and. energy.le. 3E-2) THEN
                compute_yield = eval_polynomial(coefficients_4,energy)

        ELSE IF(energy.gt. 3E-2 .and. energy.le. 5E-2) THEN 
                compute_yield = eval_polynomial(coefficients_5,energy)

    END IF 

END FUNCTION compute_yield 

   !---------------------------------------------------------------------------
   !> @author
   !> Salomon Guinchard EPFL/SPC
   !
   ! DESCRIPTION
   !> Gives random values distributed
   !> following a Poisson distrib. of parameter
   !> lambda = yield(E) for incomin ion energy E 
   !> and making use of the random_distr module (random_distr.f90)
   !--------------------------------------------------------------------------
INTEGER FUNCTION gen_elec(lambda, kmax)

    
    REAL(KIND = db) :: lambda !< Lambda parameter for Poisson distribution 
    REAL(KIND = db) :: rand   !< random number unif. generated in [0,1]
    INTEGER :: kmax           !< max number possible from Poisson
    REAL(KIND = db), DIMENSION(kmax) :: vect, SumPart, diff !< terms, partial sums for CDF and distance of rand# to CDF
    REAL(KIND = db) :: CumulPoisson !< Flag to ensure CDF ~ 1
    INTEGER :: i  
    REAL(KIND = db), DIMENSION(kmax) :: alea 

    DO i = 1,kmax
       vect(i)    = exp(-lambda)*lambda**(i-1)/factorial_fun(i-1);
       SumPart(i) = sum(vect(1:i));
       alea(i)    = 1.0 
    END DO
 
    CumulPoisson = sum(vect)
    alea = random_number*alea
    diff = abs(alea-SumPart)
    ! REMAINS TO IMPLEMENT : FIND MIN OF DIFF AND ASSOCIATE CORRECT VALUE
    ! THEN gen_elec = this value = # of generate electrons  




   ! Below: see Intel oneAPI Math Kernel Library - Fortran
   ! to optimise running speed when generating random numbers
   ! TO DO : change if enough time 
   !----------------------------------------------------------- 
   !
   ! INTEGER, INTENT(IN) :: method
   ! TYPE (VSL_STREAM_STATE), INTENT(IN) :: stream
   ! INTEGER, INTENT(IN) :: n = 1
   ! INTEGER, INTENT(IN) :: r
   ! status = virngpoisson( method, stream, n, r, lambda )
   ! gen_elec = r
   !
   !------------------------------------------------------------

END FUNCTION gen_elec






INTEGER FUNCTION factorial_fun(n)
    INTEGER :: ii 
    INTEGER :: n
    factorial_fun = 1
    IF (n .lt. 0) THEN
        RETURN
    ELSE IF (n == 0) THEN
        factorial_fun = 1
    ELSE
        DO ii=1,n
         factorial_fun = factorial_fun*ii
        END DO
   END IF
END FUNCTION factorial_fun





END MODULE iiee

