!------------------------------------------------------------------------------
! EPFL/Swiss Plasma Center
!------------------------------------------------------------------------------
!
! MODULE: IIEE 
!
!> @author
!> S. Guinchard - EPFL/SPC
!
!> Last modif.
!> 11/28 2022
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
       !#include "mkl_vsl.f90" ! for random # generators using MKL intel library

        IMPLICIT NONE

CONTAINS
        
   !---------------------------------------------------------------------------
   ! SUBROUTINE ion_induced(pion, losthole, pelec, nblostparts)  
   !
   !
   ! DESCRIPTION
   !> function to determine the number of electrons
   !> to add to a given species as a function of th number
   !> of lost ions
   !
   !--------------------------------------------------------------------------
SUBROUTINE ion_induced(pion, losthole, pelec, nblostparts)
    USE geometry 

    TYPE(particles), INTENT(INOUT):: pion, pelec !< ion and electrons parts
    REAL(KIND = db), DIMENSION(3) :: last_pos    !< last position for lost ion (revert push)
    REAL(KIND = db), DIMENSION(3) :: normal_dir  !< normal direction vector (normalised)
    INTEGER, DIMENSION(pion%Nploc):: losthole    !< indices of lost ions
    INTEGER ::i,j, nblostparts, Nploc, Nploc_old, Nploc_init!< loop indices and #lost particles
    INTEGER :: parts_size_increase, nbadded 
 
    INTEGER :: neuttype_id         !< neutral gas type_id
    INTEGER :: material_id         !< electrode material type_id 
    INTEGER :: gen_el, kmax        !< # of electrons generated, max# possibly gen. elec. 
 
    REAL(KIND = db) :: lambda             !< Poisson param. to gen elec. (yield)
    REAL(KIND = db) :: kappa, theta, Emax !< gamma distribution parameters    
    REAL(KIND = db) :: Ekin, Eem          !< kinetic energy of lost particles (yield param) and of emitted electrons
   
    kmax  = 12    !> Max num. elec. to be generated (Poisson)
    kappa = 4.0   !> kappa param. (Gamma)
    theta = 0.5   !> theta param. (Gamma) 
    Emax  = 25    !> Max value for el. (Gamma)
    
    write(*,*) 'Entering ion_induced' 
    neuttype_id = pion%neuttype_id !> temporarily stored in particle type 
    material_id = pion%material_id !> temporarily stored in particle type

    IF(pelec%Nploc + 2*nblostparts .gt. size(pelec%Z,1)) THEN
      parts_size_increase=Max(floor(0.1*size(pelec%Z,1)),2*nblostparts)
      CALL change_parts_allocation(pelec, parts_size_increase) 
    END IF
!    write(*,*) 'Error before enetering loop'
    Nploc_init = pelec%Nploc

    DO  i=1,nblostparts
      Ekin    = compute_Ekin( (/pion%UR(losthole(i)), pion%UTHET(losthole(i)), pion%UZ(losthole(i))/), pion)
      lambda  = compute_yield(Ekin, neuttype_id, material_id)
      nbadded = gen_elec(lambda, kmax)              
      Nploc_old   = pelec%Nploc
      pelec%Nploc = pelec%Nploc + nbadded
      Nploc = pelec%Nploc 
      last_pos = revert_push(pion, losthole(i))
      pelec%nbadded = pelec%nbadded+nbadded
      normal_dir = find_normal(last_pos)

      DO j=1,nbadded
        pelec%R(Nploc_old+j)     = last_pos(1)
        pelec%THET(Nploc_old+j)  = last_pos(2)
        pelec%Z(Nploc_old+j)     = last_pos(3)
        pelec%newindex = pelec%newindex + 1
        pelec%partindex(Nploc_old+j) = pelec%newindex

        IF(pelec%zero_vel == .false.) THEN
            Eem = gen_E_gamma(kappa, theta, Emax) !> generate an energy value following gamma distribution 
            pelec%UR(Nploc_old+j)    = compute_Vnorm(Eem, pelec)* normal_dir(1) !> Vr 
            pelec%UZ(Nploc_old+j)    = compute_Vnorm(Eem, pelec)* normal_dir(2) !> Vthet
            pelec%UTHET(Nploc_old+j) = compute_Vnorm(Eem, pelec)* normal_dir(3) !> Vz
        ELSE 
            pelec%UR(Nploc_old+j)    = 0.0
            pelec%UZ(Nploc_old+j)    = 0.0
            pelec%UTHET(Nploc_old+j) = 0.0
        END IF 
      END DO 
    END DO  
 !   write(*,*) 'Error  when calling geomweight'
 !   write(*,*) Nploc_init, pelec%Nploc
   if (Nploc_init-pelec%Nploc .ge. 1) then 
           call geom_weight(pelec%Z(Nploc_init+1:pelec%Nploc), pelec%R(Nploc_init+1:pelec%Nploc), pelec%geomweight(Nploc_init+1:pelec%Nploc, :)) 
   end if 
write(*,*) 'Exitting ion_induced'

END SUBROUTINE ion_induced



   !---------------------------------------------------------------------------
   ! FUNCTION compute_Ekin(velocity, p)
   !
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
   ! FUNCTION compute_Vnorm(Ekin,p)
   ! 
   !
   ! DESCRIPTION
   !> Computes the normal velocity of an incident electron emitted
   !> with energy Ekin
   !--------------------------------------------------------------------------
FUNCTION compute_Vnorm(Ekin, p) RESULT(Vnorm)
    REAL(KIND = db) :: Ekin, Vnorm !> Ekin of emitted electron, Normal. corres. veloc. 
    TYPE(particles) :: p           !> electrons      
    Vnorm = sqrt(2/p%m * Ekin * elchar) / vlight !> * elchar to get the enery in J and Vnorm in m/s
END FUNCTION compute_Vnorm



   !---------------------------------------------------------------------------
   ! FUNCTION fin_normal(last_position)
   !
   !
   ! DESCRIPTION
   !> Computes the normal velocity of an incident electron emitted
   !> with energy Ekin
   !--------------------------------------------------------------------------
FUNCTION find_normal(last_position) RESULT(normal_dir)
    USE geometry 
    REAL(KIND = db), DIMENSION(3) :: last_position !> Last pos. to eval. geom. weight at
    REAL(KIND = db), DIMENSION(3) :: normal_dir    !> Normal direction vector (Result)
    REAL(KIND = db), DIMENSION(3) :: weight        !> Geom. weight at last pos.
    REAL(KIND = db) :: norm                        !> To normalise normal vect.

    call geom_weight(last_position(3), last_position(1), weight)
    norm = sqrt(weight(2)**2 + weight(3)**2)
    normal_dir(1) = 1/norm * weight(3)  !> Normal along r
    normal_dir(2) = 0.0                 !> Normal along theta
    normal_dir(3) = 1/norm * weight(2)  !> Normal along z
END FUNCTION find_normal



   !---------------------------------------------------------------------------
   ! FUNCTION revert_push(pion, partid)
   !
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
    
    revert_push(1)  = pion%R(partid) - pion%UR(partid)*dt 
    revert_push(2)  = pion%THET(partid) -1/pion%R(partid)* pion%UTHET(partid)*dt
    revert_push(3)  = pion%Z(partid) -pion%UZ(partid)*dt 

END FUNCTION revert_push 



   !---------------------------------------------------------------------------
   ! FUNCTION eval_polynomial(coefficients, valeur)
   !
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
   ! FUNCTION compute_yield(energy, neuttype_id, material_id)
   !
   !
   ! DESCRIPTION
   !> Gives the theoretical value for the electron yield 
   !> as a function of the energy of the incident ion and
   !> the type of neutral gas  
   !
   !--------------------------------------------------------------------------
REAL(KIND = db) FUNCTION compute_yield(energy, neuttype_id, material_id) !add material id asap
    REAL(KIND = db) :: energy
    INTEGER :: neuttype_id, material_id
    REAL(KIND = db) :: Lambda_exp
    
    Lambda_exp = 1E-3 
    SELECT CASE(material_id)
        CASE(1) !304 stainless steel
                IF(energy.le. 1E-3 ) THEN 
                        SELECT CASE(neuttype_id)
                          CASE(1)
                                compute_yield = eval_polynomial(coefficients_1H_SS, energy)
                          CASE(2)
                                compute_yield = eval_polynomial(coefficients_1He_SS, energy)
                          CASE(3)
                                compute_yield = eval_polynomial(coefficients_1Ne_SS, energy)
                         CASE DEFAULT
                                compute_yield = eval_polynomial(coefficients_1H_SS, energy)
       
                        END SELECT
                ELSE IF(energy.gt. 1E-3 .and. energy.le. 1E-2) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_2_SS,energy)

                ELSE IF(energy.gt. 1E-2  .and. energy.le. 2E-2 ) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_3_SS,energy)

                ELSE IF(energy.gt. 2E-2 .and. energy.le. 3E-2) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_4_SS,energy)

                ELSE IF(energy.gt. 3E-2 .and. energy.le. 5E-2) THEN 
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_5_SS,energy)

                END IF

       CASE(2) ! Copper 
                IF(energy.le. 1E-3 ) THEN
                        SELECT CASE(neuttype_id)
                          CASE(1)
                                compute_yield = eval_polynomial(coefficients_1H_Cu, energy)
                          CASE(2)
                                compute_yield = eval_polynomial(coefficients_1He_Cu, energy)
                          CASE(3)
                                compute_yield = eval_polynomial(coefficients_1Ne_Cu, energy)
                         CASE DEFAULT
                                compute_yield = eval_polynomial(coefficients_1H_Cu, energy)

                        END SELECT
                ELSE IF(energy.gt. 1E-3 .and. energy.le. 1E-2) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_2_Cu,energy)

                ELSE IF(energy.gt. 1E-2  .and. energy.le. 2E-2 ) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_3_Cu,energy)

                ELSE IF(energy.gt. 2E-2 .and. energy.le. 3E-2) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_4_Cu,energy)

                ELSE IF(energy.gt. 3E-2 .and. energy.le. 5E-2) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_5_Cu,energy)

                END IF

       CASE(3) ! Alumium
                IF(energy.le. 1E-3 ) THEN
                        SELECT CASE(neuttype_id)
                          CASE(1)
                                compute_yield = eval_polynomial(coefficients_1H_Al, energy)
                          CASE(2)
                                compute_yield = eval_polynomial(coefficients_1He_Al, energy)
                          CASE(3)
                                compute_yield = eval_polynomial(coefficients_1Ne_Al, energy)
                         CASE DEFAULT
                                compute_yield = eval_polynomial(coefficients_1H_Al, energy)

                        END SELECT
                ELSE IF(energy.gt. 1E-3 .and. energy.le. 1E-2) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_2_Al,energy)

                ELSE IF(energy.gt. 1E-2  .and. energy.le. 2E-2 ) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_3_Al,energy)

                ELSE IF(energy.gt. 2E-2 .and. energy.le. 3E-2) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_4_Al,energy)

                ELSE IF(energy.gt. 3E-2 .and. energy.le. 5E-2) THEN
                        compute_yield = Lambda_exp * eval_polynomial(coefficients_5_Al,energy)

                END IF


    END SELECT 
END FUNCTION compute_yield 



   !---------------------------------------------------------------------------
   ! FUNCTION gen_elec(lambda, kmax)
   !
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
    INTEGER         :: kmax         !< max number possible from Poisson
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
   ! FUNCTION gen_E_gamma(kappa, theta, Emax)
   !
   !
   ! DESCRIPTION
   !> Gives random values distributed
   !> following a Gamma distrib. of parameters (kappa, theta)
   !> in [0, Emax] eV and peaked at E=2eV
   !--------------------------------------------------------------------------

FUNCTION gen_E_gamma(kappa, theta, Emax) RESULT(E_el)
 
    USE random 
    USE incomplete_gamma
        
    REAL(KIND = db) :: E_el, Emin, Emax
    REAL(KIND = db) :: kappa, theta !> parameters to shape Gamma_distr  
    REAL(KIND = db) :: nb_alea(1:1)
    REAL(KIND = db), DIMENSION(1000) :: Einc, CDF, Diff
    INTEGER :: ii !> loop index
    INTEGER :: ifault 
    INTEGER :: posid 

    Emin = 0.00 
    Emax = 25.0
    kappa = 4.0
    theta = 0.5
    
    DO ii= 1,size(Einc) 
       Einc(ii) = ( Emin + (ii-1)*(Emax-Emin)/(size(Einc)-1) )
       CDF(ii)  = gamain(Einc(ii)/theta, kappa, ifault)
    END DO 

    !> Generate Gamma distrib. int. (see Matlab. code for convg.)
    call random_array(nb_alea,1,ran_index(1),ran_array(:,1))
    Diff  = abs(nb_alea(1) - CDF)
    posid = minloc(Diff, dim = 1)

        
    E_el  = Einc(posid)
END FUNCTION gen_E_gamma

   !---------------------------------------------------------------------------
   ! FUNCTION factorial_fun(n)
   !
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

