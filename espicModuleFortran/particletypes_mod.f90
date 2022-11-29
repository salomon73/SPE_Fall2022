!------------------------------------------------------------------------------
! EPFL/Swiss Plasma Center
!------------------------------------------------------------------------------
!
! MODULE: particletypes
!
!> @author
!> Guillaume Le Bars EPFL/SPC
!> Patryk Kaminski   EPFL/SPC
!> Trach Minh Tran   EPFL/SPC
!
! DESCRIPTION:
!> Module responsible for defining the particle types and defining some subroutines to change their size,
!> initialize them or delete them
!------------------------------------------------------------------------------

MODULE particletypes
    USE constants
    !
      IMPLICIT NONE
    
      !> Stores the particles properties for the run.
      TYPE particles
        INTEGER :: Nploc                                     !< Local number of simulated particles
        INTEGER :: Nptot                                     !< Total number of simulated particles
        INTEGER :: Newindex                                  !< Stores the higher partindex for the creation of new particles
        REAL(kind=db) :: m                                   !< Particle mass
        REAL(kind=db) :: q                                   !< Particle charge
        REAL(kind=db) :: weight                              !< Number of particles represented by one macro-particle
        REAL(kind=db) :: qmRatio                             !< Charge over mass ratio
        REAL(kind=db) :: nudcol(3)                           !< Effective momentum drag frequency 
        REAL(kind=db) :: H0
        REAL(kind=db) :: P0
        REAL(kind=db) :: temperature
        LOGICAL :: Davidson=.false.
        LOGICAL :: is_test= .false.                          !< determines if particle is saved on ittracer
        LOGICAL :: is_field= .true.                          !< determines if particle contributes to Poisson solver
        LOGICAL :: calc_moments=.false.
        INTEGER, allocatable               :: nblost(:)      !< number of particles lost in domain boundaries at current timestep
        INTEGER                            :: nbadded        !< number of particles added by source since last gather
        INTEGER, DIMENSION(2)              :: nbcolls        !< number of particles collisions with neutrals ionisation,   elastic)
        INTEGER, DIMENSION(:), ALLOCATABLE :: Rindex         !< Index in the electric potential grid for the R direction
        INTEGER, DIMENSION(:), ALLOCATABLE :: Zindex         !< Index in the electric potential grid for the Z direction
        INTEGER, DIMENSION(:), ALLOCATABLE :: partindex      !< Index of the particle to be able to follow it when it goes from one MPI host to the other
        INTEGER                            :: iiee_id=-1     !< Index defining whether or not ion induced ee are considered
        INTEGER                            :: neuttype_id=1  !< Index defining which type of neutral gas is used to produce the ions
        INTEGER                            :: material_id=1  !< Index defining the type of material for the electrodes (1=304SS)
        LOGICAL                            :: zero_vel=.true.!< Defines wether or not the electrons are gen. with init. vel
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: R     !< radial coordinates of the particles
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: Z     !< longitudinal coordinates of the particles
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: THET  !< azimuthal coordinates of the particles
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: BZ    !< axial radial relative distances to the left grid line
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: BR    !< radial relative distances to the bottom grid line
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: pot   !< Electric potential
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: potxt !< External electric potential
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: Er    !< Radial Electric field
        REAL(kind=db), DIMENSION(:), ALLOCATABLE :: Ez    !< Axial electric field
        REAL(kind=db), DIMENSION(:), CONTIGUOUS, POINTER:: UR          !< normalized radial velocity at the current time step
        REAL(kind=db), DIMENSION(:), CONTIGUOUS, POINTER:: URold       !< normalized radial velocity at the previous time step
        REAL(kind=db), DIMENSION(:), CONTIGUOUS, POINTER:: UTHET       !< normalized azimuthal velocity at the current time step
        REAL(kind=db), DIMENSION(:), CONTIGUOUS, POINTER:: UTHETold    !< normalized azimuthal velocity at the previous time step
        REAL(kind=db), DIMENSION(:), CONTIGUOUS, POINTER:: UZ          !< normalized axial velocity at the current time step
        REAL(kind=db), DIMENSION(:), CONTIGUOUS, POINTER:: UZold       !< normalized axial velocity at the previous time step
        REAL(kind=db), DIMENSION(:), CONTIGUOUS, POINTER:: Gamma       !< Lorentz factor at the current time step
        REAL(kind=db), DIMENSION(:), CONTIGUOUS, POINTER:: Gammaold    !< Lorentz factor at the previous time step
        Real(kind=db), Dimension(:,:),ALLOCATABLE:: geomweight !< geometric weight at the particle position
        Real(kind=db), Dimension(:,:),ALLOCATABLE:: moments !< stores the moment matrix
        LOGICAL:: collected                               !< Stores if the particles data have been collected to MPI root process during this timestep
        INTEGER, DIMENSION(:), ALLOCATABLE:: addedlist
      END TYPE particles

      !> Structure containing a single particle position and velocity used in MPI communications.
      TYPE particle
        INTEGER ::       partindex =0
        REAL(kind=db) :: R         =0
        REAL(kind=db) :: Z         =0
        REAL(kind=db) :: THET      =0
        REAL(kind=db) :: UZ        =0
        REAL(kind=db) :: UR        =0
        REAL(kind=db) :: UTHET     =0
        REAL(kind=db) :: Gamma     =0
        REAL(kind=db) :: pot       =0
      END TYPE particle
    
      TYPE linked_part
        type(particle) p
        type(linked_part), POINTER:: next=> NULL()
        type(linked_part), POINTER:: prev=> NULL()
      END TYPE linked_part
    
      TYPE linked_part_row
        INTEGER :: n = 0
        type(linked_part), POINTER:: start=>NULL()
        type(linked_part), POINTER:: end=>NULL()
      END TYPE linked_part_row
  
      CONTAINS

      !---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Allocate the memory for the particles variable storing the particles quantities.
!
!> @param[inout] p the particles variable needing to be allocated.
!> @param[in] nparts the maximum number of particles that will be stored in this variable
!---------------------------------------------------------------------------

SUBROUTINE creat_parts(p, nparts)
  TYPE(particles)       :: p
  INTEGER, INTENT(in)   :: nparts

  IF (.NOT. ALLOCATED(p%Z) ) THEN
    p%Nploc = nparts
    p%Nptot = nparts
    ALLOCATE(p%Z(nparts))
    ALLOCATE(p%R(nparts))
    ALLOCATE(p%THET(nparts))
    ALLOCATE(p%BZ(nparts))
    ALLOCATE(p%BR(nparts))
    ALLOCATE(p%UR(nparts))
    ALLOCATE(p%UZ(nparts))
    ALLOCATE(p%UTHET(nparts))
    ALLOCATE(p%URold(nparts))
    ALLOCATE(p%UZold(nparts))
    ALLOCATE(p%UTHETold(nparts))
    ALLOCATE(p%Gamma(nparts))
    ALLOCATE(p%Rindex(nparts))
    ALLOCATE(p%Zindex(nparts))
    ALLOCATE(p%partindex(nparts))
    ALLOCATE(p%pot(nparts))
    ALLOCATE(p%potxt(nparts))
    ALLOCATE(p%Er(nparts))
    ALLOCATE(p%Ez(nparts))
    ALLOCATE(p%GAMMAold(nparts))
    Allocate(p%geomweight(nparts,0:2))
    if(.not. allocated(p%nblost) )then
        allocate(p%nblost(4))
    end if
    p%newindex=0
    p%nblost=0
    p%nbadded=0
    p%partindex=-1
    p%iiee_id=-1
    p%neuttype_id=1
    p%material_id=1
    p%zero_vel=.true.
    p%URold=0
    p%UZold=0
    p%UTHETold=0
    p%rindex=0
    p%zindex=0
    p%BR=0
    p%BZ=0
    p%UR=0
    p%UZ=0
    p%UTHET=0
    p%Z=0
    p%R=0
    p%THET=0
    p%Gamma=1
    p%Er=0
    p%Ez=0
    p%pot=0
    p%potxt=0
    p%gammaold=1
    p%collected=.false.
    p%Davidson=.false.
    p%is_test=.false.
    p%is_field=.true.
    p%calc_moments=.true.
    p%m=me
    p%q=-elchar
    p%qmRatio=p%q/p%m
    p%weight=1.0_db
    p%H0=0
    p%P0=0
    p%temperature=0
    p%geomweight=0
  END IF
END SUBROUTINE creat_parts

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Copy one particle from the receive buffers to the local simulation variable parts.
!
!> @param [in] part particle parameters to copy from
!> @param [in] partsindex destination particle index in the local parts variable
!---------------------------------------------------------------------------
SUBROUTINE Insertincomingpart(p, part, partsindex)
  TYPE(particles), INTENT(INOUT):: p
  INTEGER, INTENT(in) :: partsindex
  TYPE(particle), INTENT(in) :: part
    p%partindex(partsindex) = part%partindex
    p%R(partsindex) =         part%R
    p%Z(partsindex) =         part%Z
    p%THET(partsindex) =      part%THET
    p%UZ(partsindex) =        part%UZ
    p%UR(partsindex) =        part%UR
    p%UTHET(partsindex) =     part%UTHET
    p%Gamma(partsindex) =     part%Gamma
    p%pot(partsindex) =       part%pot
!
END SUBROUTINE Insertincomingpart

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Copy one particle from the local parts variable to the send buffer.
!
!> @param [in] buffer send buffer to copy to
!> @param [in] bufferindex particle index in the send buffer
!> @param [in] partsindex origin particle index in the local parts variable
!---------------------------------------------------------------------------
SUBROUTINE Insertsentpart(p, buffer, bufferindex, partsindex)
  TYPE(particles), INTENT(INOUT):: p
  INTEGER, INTENT(in) :: bufferindex, partsindex
  TYPE(particle), DIMENSION(:), INTENT(inout) :: buffer
    buffer(bufferindex)%partindex = p%partindex(partsindex)
    buffer(bufferindex)%R         = p%R(partsindex)
    buffer(bufferindex)%Z         = p%Z(partsindex)
    buffer(bufferindex)%THET      = p%THET(partsindex)
    buffer(bufferindex)%UZ        = p%UZ(partsindex)
    buffer(bufferindex)%UR        = p%UR(partsindex)
    buffer(bufferindex)%UTHET     = p%UTHET(partsindex)
    buffer(bufferindex)%Gamma     = p%Gamma(partsindex)
    buffer(bufferindex)%pot       = p%pot(partsindex)
!
END SUBROUTINE Insertsentpart


!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief Exchange two particles in the parts variable.
!
!> @param [in] index1 index in parts of the first particle to exchange.
!> @param [in] index2 index in parts of the second particle to exchange.
!---------------------------------------------------------------------------
SUBROUTINE exchange_parts(p, index1, index2)
  TYPE(particles), INTENT(INOUT):: p
  INTEGER, INTENT(IN) :: index1, index2
  REAL(kind=db):: R, Z, THET, UR, UZ, UTHET, Gamma, geomweight(0:2),pot
  INTEGER :: Rindex, Zindex, partindex
  !! Exchange particle at index1 with particle at index2

  ! Store part at index1 in temporary value
  partindex  = p%partindex(index1)
  Gamma      = p%Gamma(index1)
  pot        = p%pot(index1)
  R          = p%R(index1)
  Z          = p%Z(index1)
  THET       = p%THET(index1)
  UR         = p%UR(index1)
  UTHET      = p%UTHET(index1)
  UZ         = p%UZ(index1)
  Rindex     = p%Rindex(index1)
  Zindex     = p%Zindex(index1)
  geomweight = p%geomweight(index1,:)

  ! Move part at index2 in part at index 1
  p%partindex(index1)    = p%partindex(index2)
  p%Gamma(index1)        = p%Gamma(index2)
  p%pot(index1)          = p%pot(index2)
  p%R(index1)            = p%R(index2)
  p%Z(index1)            = p%Z(index2)
  p%THET(index1)         = p%THET(index2)
  p%UR(index1)           = p%UR(index2)
  p%UTHET(index1)        = p%UTHET(index2)
  p%UZ(index1)           = p%UZ(index2)
  p%Rindex(index1)       = p%Rindex(index2)
  p%Zindex(index1)       = p%Zindex(index2)
  p%geomweight(index1,:) = p%geomweight(index2,:)

  ! Move temporary values from part(index1) to part(index2)
  p%partindex(index2)    = partindex
  p%Gamma(index2)        = Gamma
  p%pot(index2)          = pot
  p%R(index2)            = R
  p%Z(index2)            = Z
  p%THET(index2)         = THET
  p%UR(index2)           = UR
  p%UTHET(index2)        = UTHET
  p%UZ(index2)           = UZ
  p%Rindex(index2)       = Rindex
  p%Zindex(index2)       = Zindex
  p%geomweight(index2,:) = geomweight

  END SUBROUTINE exchange_parts

  SUBROUTINE change_parts_allocation(p, sizedifference)
    implicit none
    TYPE(particles), INTENT(INOUT):: p
    INTEGER,INTENT(IN) :: sizedifference
    CALL change_array_size_int(p%Rindex, sizedifference)
    CALL change_array_size_int(p%Zindex, sizedifference)
    CALL change_array_size_int(p%partindex, sizedifference)
    CALL change_array_size_dp(p%ER,sizedifference)
    CALL change_array_size_dp(p%EZ,sizedifference)
    CALL change_array_size_dp(p%pot,sizedifference)
    CALL change_array_size_dp(p%potxt,sizedifference)
    CALL change_array_size_dp(p%R,sizedifference)
    CALL change_array_size_dp(p%Z,sizedifference)
    CALL change_array_size_dp(p%THET,sizedifference)
    CALL change_array_size_dp(p%BR,sizedifference)
    CALL change_array_size_dp(p%BZ,sizedifference)
    CALL change_array_size_dp2(p%geomweight,sizedifference)
    CALL change_array_size_dp_ptr(p%UR,sizedifference)
    CALL change_array_size_dp_ptr(p%URold,sizedifference)
    CALL change_array_size_dp_ptr(p%UZ,sizedifference)
    CALL change_array_size_dp_ptr(p%UZold,sizedifference)
    CALL change_array_size_dp_ptr(p%UTHET,sizedifference)
    CALL change_array_size_dp_ptr(p%UTHETold,sizedifference)
    CALL change_array_size_dp_ptr(p%Gamma,sizedifference)
    CALL change_array_size_dp_ptr(p%Gammaold,sizedifference)
    p%Nploc=MIN(p%Nploc,size(p%R))
  END SUBROUTINE change_parts_allocation

  SUBROUTINE change_array_size_dp(arr, sizedifference)
    implicit none
    REAL(kind=db), ALLOCATABLE, INTENT(INOUT):: arr(:)
    INTEGER, INTENT(IN):: sizedifference
    REAL(kind=db), ALLOCATABLE:: temp(:)
    INTEGER:: current_size, new_size
    if(allocated(arr)) THEN
      current_size=size(arr)
      new_size=current_size+sizedifference
      ALLOCATE(temp(new_size))
      temp(1:min(current_size,new_size))=arr(1:min(current_size,new_size))
      DEALLOCATE(arr)
      CALL move_alloc(temp, arr)
    END IF
  END SUBROUTINE change_array_size_dp

  SUBROUTINE change_array_size_dp2(arr, sizedifference)
    implicit none
    REAL(kind=db), ALLOCATABLE, INTENT(INOUT):: arr(:,:)
    INTEGER, INTENT(IN):: sizedifference
    REAL(kind=db), ALLOCATABLE:: temp(:,:)
    INTEGER:: current_size, new_size
    if(allocated(arr)) THEN
      current_size=size(arr,1)
      new_size=current_size+sizedifference
      ALLOCATE(temp(new_size,0:size(arr,2)-1))
      temp(1:min(current_size,new_size),:)=arr(1:min(current_size,new_size),:)
      DEALLOCATE(arr)
      CALL move_alloc(temp, arr)
    END IF
  END SUBROUTINE change_array_size_dp2

  SUBROUTINE change_array_size_dp_ptr(arr, sizedifference)
    implicit none
    REAL(kind=db), POINTER, INTENT(INOUT):: arr(:)
    INTEGER, INTENT(IN):: sizedifference
    REAL(kind=db),  CONTIGUOUS, POINTER:: temp(:)
    INTEGER:: current_size, new_size
    if(associated(arr)) THEN
      current_size=size(arr)
      new_size=current_size+sizedifference
      ALLOCATE(temp(new_size))
      temp(1:min(current_size,new_size))=arr(1:min(current_size,new_size))
      DEALLOCATE(arr)
      arr=> temp
    END IF
  END SUBROUTINE change_array_size_dp_ptr

  SUBROUTINE change_array_size_int(arr, sizedifference)
    implicit none
    INTEGER, ALLOCATABLE, INTENT(INOUT):: arr(:)
    INTEGER, INTENT(IN):: sizedifference
    INTEGER, ALLOCATABLE:: temp(:)
    INTEGER:: current_size, new_size

    if(allocated(arr)) THEN
      current_size=size(arr)
      new_size=current_size+sizedifference
      ALLOCATE(temp(new_size))
      temp(1:min(current_size,new_size))=arr(1:min(current_size,new_size))
      DEALLOCATE(arr)
      CALL move_alloc(temp,arr)
    END IF
  END SUBROUTINE change_array_size_int

  !---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Move particle with index sourceindex to particle with index destindex.
!> !WARNING! This will overwrite particle at destindex.
!
!> @param [in] sourceindex index in parts of the particle to move.
!> @param [in] destindex index in parts of the moved particle destination.
!---------------------------------------------------------------------------
  SUBROUTINE move_part(p, sourceindex, destindex)
    !! This will destroy particle at destindex
    INTEGER, INTENT(IN) :: destindex, sourceindex
    TYPE(particles), INTENT(INOUT)::p
  
    IF(sourceindex .eq. destindex) RETURN
    IF(sourceindex .le. 0 .or. destindex .le. 0) RETURN
    ! Move part at sourceindex in part at destindex
   Call copy_part(p,sourceindex,destindex,p)
  
  END SUBROUTINE move_part

  !---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Copy particle with index sourceindex in particles sourcep to particle with index destindex in particles destp.
!> !WARNING! This will overwrite particle at destp(destindex).
!
!> @param [inout] sourcep Structure of source particles.
!> @param [in] sourceindex index in parts of the particle to move.
!> @param [in] destindex index in parts of the moved particle destination.
!> @param [inout] destp Structure of source particles.
!---------------------------------------------------------------------------
    SUBROUTINE copy_part(sourcep, sourceindex, destindex, destp)
      !! This will destroy particle at destindex
      INTEGER, INTENT(IN) :: destindex, sourceindex
      TYPE(particles), INTENT(IN)::sourcep
      TYPE(particles), INTENT(INOUT)::destp
    
      IF(sourceindex .le. 0 .or. destindex .le. 0) RETURN
      IF( destindex .gt. size(destp%R,1)) RETURN
      ! Move part at sourceindex in part at destindex
      destp%partindex(destindex)    = sourcep%partindex(sourceindex)
      destp%Gamma(destindex)        = sourcep%Gamma(sourceindex)
      destp%Gammaold(destindex)     = sourcep%Gammaold(sourceindex)
      destp%R(destindex)            = sourcep%R(sourceindex)
      destp%Z(destindex)            = sourcep%Z(sourceindex)
      destp%THET(destindex)         = sourcep%THET(sourceindex)
      destp%UR(destindex)           = sourcep%UR(sourceindex)
      destp%UTHET(destindex)        = sourcep%UTHET(sourceindex)
      destp%UZ(destindex)           = sourcep%UZ(sourceindex)
      destp%URold(destindex)        = sourcep%URold(sourceindex)
      destp%UTHETold(destindex)     = sourcep%UTHETold(sourceindex)
      destp%UZold(destindex)        = sourcep%UZold(sourceindex)
      destp%Rindex(destindex)       = sourcep%Rindex(sourceindex)
      destp%Zindex(destindex)       = sourcep%Zindex(sourceindex)
      destp%geomweight(destindex,:) = sourcep%geomweight(sourceindex,:)
      destp%pot(destindex)          = sourcep%pot(sourceindex)
      destp%potxt(destindex)        = sourcep%potxt(sourceindex)
    
    END SUBROUTINE copy_part
  !________________________________________________________________________________

    SUBROUTINE destroy_parts(p)
      TYPE(particles)       :: p
      p%Nploc=0
      IF(ALLOCATED(p%Z)) DEALLOCATE(p%Z)
      IF(ALLOCATED(p%R)) DEALLOCATE(p%R)
      IF(ALLOCATED(p%THET)) DEALLOCATE(p%THET)
      IF(ALLOCATED(p%BZ)) DEALLOCATE(p%BZ)
      IF(ALLOCATED(p%BR)) DEALLOCATE(p%BR)
      IF(ASSOCIATED(p%UR)) DEALLOCATE(p%UR)
      IF(Associated(p%URold)) DEALLOCATE(p%URold)
      IF(Associated(p%UZ)) DEALLOCATE(p%UZ)
      IF(Associated(p%UZold)) DEALLOCATE(p%UZold)
      IF(Associated(p%UTHET)) DEALLOCATE(p%UTHET)
      IF(Associated(p%UTHETold)) DEALLOCATE(p%UTHETold)
      IF(Associated(p%Gamma)) DEALLOCATE(p%Gamma)
      IF(Associated(p%Gammaold)) DEALLOCATE(p%Gammaold)
      IF(ALLOCATED(p%Rindex)) DEALLOCATE(p%Rindex)
      IF(ALLOCATED(p%Zindex)) DEALLOCATE(p%Zindex)
      IF(ALLOCATED(p%partindex)) DEALLOCATE(p%partindex)
      if(allocated(p%geomweight)) Deallocate(p%geomweight)
      if(allocated(p%moments)) Deallocate(p%moments)
    END SUBROUTINE
  !________________________________________________________________________________
    SUBROUTINE clean_beam(partslist)
  !
      INTEGER:: i
      type(particles):: partslist(:)

      Do i=1,size(partslist,1)
        CALL destroy_parts(partslist(i))
      END DO
  !
    END SUBROUTINE clean_beam
  !________________________________________________________________________________
    SUBROUTINE swappointer( pointer1, pointer2)
      REAL(kind=db), DIMENSION(:), POINTER, INTENT(inout):: pointer1, pointer2
      REAL(kind=db), DIMENSION(:), POINTER:: temppointer
      temppointer=>pointer1
      pointer1=>pointer2
      pointer2=>temppointer
    END SUBROUTINE swappointer


!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Deallocate recursively a linked_paticle linked list 
!
!> @param [in] l_p linked_part particle to be dallocated.
!---------------------------------------------------------------------------
  RECURSIVE SUBROUTINE destroy_linked_parts(l_p)
    TYPE(linked_part), POINTER :: l_p

    IF(associated(l_p%next)) call destroy_linked_parts(l_p%next)
    deallocate(l_p)
  END subroutine destroy_linked_parts


!--------------------------------------------------------------------------
!> @author 
!> S.Guinchard  EPFL/SPC
!> Last modified on: 11/15/2022
! 
!DESCRIPTION 
!> Function giving particle energy for a given partindex
!
!-------------------------------------------------------------------------

REAL(KIND=db) FUNCTION eKin_part(p, partind)
    TYPE(particles), INTENT(INOUT):: p
    INTEGER :: partind

        eKin_part = 0.5* p%m * (p%UR(partind)**2 + p%UTHET(partind)**2 + p%UZ(partind)**2 )

END FUNCTION eKin_part



END MODULE particletypes
