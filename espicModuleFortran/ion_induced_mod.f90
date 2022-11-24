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
        USE materials
        !USE factorial  !    ASK GUILLAUME HOW TO COMBINE MODULES 
        !#include "mkl_vsl.f90"

        IMPLICIT NONE

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
    INTEGER ::i,j, nblostparts, Nploc, Nploc_old !< loop indices and #lost particles
    INTEGER :: parts_size_increase, nbadded 
    INTEGER :: neuttype_id         !< neutral gas type_id 
    INTEGER :: gen_el, kmax        !< # of electrons generated, max# possibly gen. elec. 
    REAL(KIND = db) :: lambda      !< Poisson param. to gen elec. (yield)       
    REAL(KIND = db) :: Ekin        !< kinetic energy of lost particles (yield param)
    REAL(KIND = db) :: dr          !< dr displacement for created parts
                                   !< will be changed to normal disp       

 
    !nbadded=nblostparts !< # particles to add to electron species
    dr   = 1E-3/rnorm
    kmax = 12
    neuttype_id = pion%neuttype_id !> temporarily stored in particle type 


    IF(pelec%Nploc + 2*nblostparts .gt. size(pelec%Z,1)) THEN
      parts_size_increase=Max(floor(0.1*size(pelec%Z,1)),2*nblostparts)
      CALL change_parts_allocation(pelec, parts_size_increase) 
    END IF

    DO  i=1,nblostparts
      Ekin    = compute_Ekin( (/pion%UR(losthole(i)), pion%UTHET(losthole(i)), pion%UZ(losthole(i))/), pion)
      lambda  = compute_yield(Ekin, neuttype_id)
      nbadded = gen_elec(lambda, kmax)  
      Nploc_old   = pelec%Nploc 
      pelec%Nploc = pelec%Nploc + nbadded
      Nploc = pelec%Nploc 
      last_pos = revert_push(pion, losthole(i))
      pelec%nbadded = pelec%nbadded+nbadded
      DO j=1,nbadded
        pelec%R(Nploc_old+j)     = last_pos(1)
        pelec%THET(Nploc_old+j)  = last_pos(2)
        pelec%Z(Nploc_old+j)     = last_pos(3)
        pelec%UR(Nploc_old+j)    = 0
        pelec%UZ(Nploc_old+j)    = 0
        pelec%UTHET(Nploc_old+j) = 0
      END DO 
    END DO    
END SUBROUTINE ion_induced



   !---------------------------------------------------------------------------
   !> @author
   !> Salomon Guinchard EPFL/SPC
   !
   ! DESCRIPTION
   !> Computes the kinetic energy of a particle given its 3-vel. components
   !--------------------------------------------------------------------------
FUNCTION  compute_Ekin(velocity, p) RESULT(Ekin)

    TYPE(particles), INTENT(INOUT):: p
    REAL(KIND = db), DIMENSION(3) :: velocity
    REAL(KIND = db) :: Ekin
  
    Ekin = 5E-7 * p%m * vlight**2 /elchar * (velocity(1)**2 + velocity(2)**2 + velocity(3)**2)

END FUNCTION compute_Ekin



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
    INTEGER :: partid                     !> id of particle to reverse position               
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
REAL(KIND = db) FUNCTION eval_polynomial(coefficients, valeur)
    REAL(KIND = db), DIMENSION(:) :: coefficients !< polynomial (e.g fitted yield) coeffs
    REAL(KIND = db) :: valeur                     !< point where to evaluate polyn
    INTEGER :: ii

    eval_polynomial = 0
    DO ii=1, size(coefficients)
      eval_polynomial = eval_polynomial+coefficients(ii)*valeur**(ii-1)
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
REAL(KIND = db) FUNCTION compute_yield(energy, neuttype_id) !add material id asap
    REAL(KIND = db) :: energy
    INTEGER :: neuttype_id, material_id
    REAL(KIND = db) :: Lambda_exp
    
    Lambda_exp = 1E-3 
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
                compute_yield = Lambda_exp * eval_polynomial(coefficients_2,energy)

        ELSE IF(energy.gt. 1E-2  .and. energy.le. 2E-2 ) THEN
                compute_yield = Lambda_exp * eval_polynomial(coefficients_3,energy)

        ELSE IF(energy.gt. 2E-2 .and. energy.le. 3E-2) THEN
                compute_yield = Lambda_exp * eval_polynomial(coefficients_4,energy)

        ELSE IF(energy.gt. 3E-2 .and. energy.le. 5E-2) THEN 
                compute_yield = Lambda_exp * eval_polynomial(coefficients_5,energy)

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
   !--------------------------------------------------------------------------
INTEGER FUNCTION gen_elec(lambda, kmax)

    USE random  
    REAL(KIND = db) :: lambda       !< Lambda parameter for Poisson distribution 
    REAL(KIND = db) :: nb_alea(1:1) !< random number unif. generated in [0,1]
    INTEGER :: kmax           !< max number possible from Poisson
    REAL(KIND = db) :: CumulPoisson  !< Flag to ensure CDF ~ 1
    INTEGER :: i, ii                 !< loop indices    
    REAL(KIND = db), DIMENSION(kmax) :: vect, SumPart !< terms, partial sums for CDF
       
    !> Compute probabilities for each int. value and CDF values
    DO i = 1,kmax
       vect(i)    = exp(-lambda)*lambda**(i-1)/factorial_fun(i-1);
       SumPart(i) = sum(vect(1:i));
    END DO
    
    !> CDF expected to sum to 1 
    CumulPoisson = sum(vect)

   ! call random_number(nb_alea)
   ! DO ii = 1,size(SumPart)-1
   !     IF (nb_alea .lt. SumPart(1)) THEN
   !             gen_elec = 0
   !     ELSE IF ((SumPart(ii).le.nb_alea) .and. (nb_alea .lt. SumPart(ii+1))) THEN
   !             gen_elec = ii
   !     END IF
   ! END DO



    !> Generate poisson distrib. int. (see Matlab. code for convg.)
    call random_array(nb_alea,1,ran_index(1),ran_array(:,1))  
    DO ii = 1,size(SumPart)-1
       IF (nb_alea(1) .lt. SumPart(1)) THEN
                gen_elec = 0
        ELSE IF ((SumPart(ii).le.nb_alea(1)) .and. (nb_alea(1) .lt. SumPart(ii+1))) THEN
                gen_elec = ii  
        END IF  
     END DO 
    

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



   !---------------------------------------------------------------------------
   !> @author
   !> Salomon Guinchard EPFL/SPC
   !
   ! DESCRIPTION
   !> Gives the factorial of an integer
   !--------------------------------------------------------------------------
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

