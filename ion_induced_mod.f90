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
        pelec%R(Nploc)    = pion%R(losthole(j)) + dr
        pelec%Z(Nploc)    = pion%Z(losthole(j)) 
        pelec%THET(Nploc) = pion%THET(losthole(j))

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
SUBROUTINE invert_push(pion)
    TYPE(particles), INTENT(INOUT):: pion


END SUBROUTINE invert_push 


END MODULE iiee
