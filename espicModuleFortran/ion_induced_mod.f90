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
        IMPLICIT NONE
        REAL(KIND = db) :: yield  !< net electronic yield for incident ion
                                  !< yield is a module variable
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
    revert_push(2)  = pion%THET(partid) - pion%UTHET(partid)*dt
    revert_push(3)  = pion%Z(partid) -pion%UZ(partid)*dt 

    ! BELOW WE TRY TO REVERSE THE ANGLE 
    ! REMAINS TO BE DONE 


END FUNCTION revert_push 


END MODULE iiee

