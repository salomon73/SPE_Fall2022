MODULE beam
!------------------------------------------------------------------------------
! EPFL/Swiss Plasma Center
!------------------------------------------------------------------------------
!
! MODULE: beam
!
!> @author
!> Guillaume Le Bars EPFL/SPC
!> Patryk Kaminski   EPFL/SPC
!> Trach Minh Tran   EPFL/SPC
!
! DESCRIPTION:
!> Module responsible for loading, advancing and computing the necessary diagnostics for the simulated particles.
!------------------------------------------------------------------------------
!
  USE constants
  use mpi
  USE mpihelper
  USE basic, ONLY: mpirank, mpisize
  USE distrib
  USE particletypes
  USE weighttypes
  

  IMPLICIT NONE

!
  !TYPE(particles) :: parts  !< Storage for all the particles
  !SAVE :: parts
  TYPE(particles), DIMENSION(:), ALLOCATABLE, SAVE :: partslist
  
!   Diagnostics (scalars)
  REAL(kind=db) :: ekin=0  !< Total kinetic energy (J)
  REAL(kind=db) :: epot=0  !< Total potential energy (J)
  REAL(kind=db) :: etot=0  !< Current total energy (J)
  REAL(kind=db) :: etot0=0 !< Initial total energy (J)
  REAL(kind=db) :: loc_etot0=0 !< theoretical local total energy (J)
  REAL(kind=db) :: Energies(4) !< (1) kinetic energy, (2) potential energy, (3) total energy and (4) gained/lossed energy due to gain or loss of particles (J) 
!
 INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: Nplocs_all !< Array containing the local numbers of particles in each MPI process
 
 INTERFACE add_created_part
    MODULE PROCEDURE add_linked_created_part, add_list_created_part
 END INTERFACE add_created_part
!

 abstract interface
 subroutine rloader(nbase,y,rminus,rplus)
   USE constants
   REAL(kind=db), INTENT(out) :: y(:)
   INTEGER, INTENT(in)        :: nbase
   REAL(kind=db), INTENT(in)  :: rplus, rminus
 end subroutine
  REAL(kind=db) FUNCTION gamma(UZ, UR, UTHET)
      USE constants
      REAL(kind=db), INTENT(IN):: UR,UZ,UTHET
  end FUNCTION
end interface

CONTAINS

!---------------------------------------------------------------------------
   !> @author
   !> Guillaume Le Bars EPFL/SPC
   !
   ! DESCRIPTION:
   !> @brief Loads the particles at the beginning of the simulation and create the parts variable if necessary
   !---------------------------------------------------------------------------
SUBROUTINE load_parts
    USE basic, ONLY: nplasma, mpirank, ierr, distribtype, mpisize, nlclassical, nbspecies, Zbounds, partfile
    use mpi

    INTEGER:: i
    REAL(kind=db), DIMENSION(:), ALLOCATABLE :: VZ, VR, VTHET

    
    ALLOCATE(VZ(nplasma), VR(nplasma), VTHET(nplasma))


    ! Select case to define the type of distribution
    SELECT CASE(distribtype)
      CASE(1) ! Gaussian distribution in V, uniform in Z and 1/R in R
        CALL loaduniformRZ(partslist(1), VR, VZ, VTHET)
      CASE(2) !Stable distribution from Davidson 4.95 p.119
        CALL loadDavidson(partslist(1), VR, VZ, VTHET, lodunir)
      CASE(3) !Stable distribution from Davidson 4.95 p.119 but with constant distribution in R 
        CALL loadDavidson(partslist(1), VR, VZ, VTHET, lodinvr)
      CASE(4) !Stable distribution from Davidson 4.95 p.119 but with gaussian distribution in R
        CALL loadDavidson(partslist(1), VR, VZ, VTHET, lodgausr)
      CASE(5) !Stable distribution from Davidson 4.95 p.119 with gaussian in V computed from v_th given by temp
        CALL loadDavidson(partslist(1), VR, VZ, VTHET, lodunir)
      CASE(6) ! Uniform distribution in R and Z and Gaussian distribution in V with Vz<V_perp to satisfy magnetic mirror trapping
        CALL loaduniformRZ(partslist(1), VR, VZ, VTHET)
        VZ = VZ/50
      CASE(7) ! Distribution defined in separate input file
        CALL read_part_file(partslist(1), partfile(1), VR, VZ, VTHET) 
        nplasma=partslist(1)%Nptot
      CASE DEFAULT
        IF (mpirank .eq. 0) WRITE(*,*) "Unknown type of distribution:", distribtype
        CALL MPI_Abort(MPI_COMM_WORLD, -1, ierr)

    END SELECT
    Do i=1,nplasma
      partslist(1)%partindex(i)=i
    END DO
    partslist(1)%newindex=nplasma

    IF(nlclassical) THEN
      partslist(1)%Gamma(1:nplasma)=1.0_db
    ELSE
      partslist(1)%Gamma(1:nplasma)=sqrt(1/(1-VR**2-VZ**2-VTHET**2))
    END IF
    ! Normalization of the velocities
    partslist(1)%UR(1:nplasma)=partslist(1)%Gamma(1:nplasma)*VR
    partslist(1)%UZ(1:nplasma)=partslist(1)%Gamma(1:nplasma)*VZ
    partslist(1)%UTHET(1:nplasma)=partslist(1)%Gamma(1:nplasma)*VTHET
    DEALLOCATE(VZ, VR, VTHET)
    
    CALL boundary_loss(partslist(1))
    CALL localisation(partslist(1))
    partslist(1)%Nptot=partslist(1)%nploc

    
    DO i=2,nbspecies
     call load_part_file(partslist(i),partfile(i))
    END DO

    partslist(1)%calc_moments=.true.
    partslist(1)%is_field=.true.
END SUBROUTINE load_parts

SUBROUTINE load_part_file(p,partfilename)
  use mpi
  USE basic, ONLY: mpisize, nlclassical, Zbounds
  type(particles) :: p
  CHARACTER(len=*):: partfilename
  Real(kind=db), ALLOCATABLE:: VR(:), VZ(:), VTHET(:)
  INTEGER j
  ! Read the actual file and load the position in p
  CALL read_part_file(p, partfilename, VR, VZ, VTHET)
  Do j=1,p%Nploc
    p%partindex(j)=j
  END DO

  IF(nlclassical) THEN
    p%Gamma=1.0_db
  ELSE
    p%Gamma(1:p%Nptot)=sqrt(1/(1-VR**2-VZ**2-VTHET**2))
  END IF
  ! Normalization of the velocities
  p%UR(1:p%Nptot)   = p%Gamma(1:p%Nptot)*VR(1:p%Nptot)
  p%UZ(1:p%Nptot)   = p%Gamma(1:p%Nptot)*VZ(1:p%Nptot)
  p%UTHET(1:p%Nptot)= p%Gamma(1:p%Nptot)*VTHET(1:p%Nptot)
  DEALLOCATE(VZ, VR, VTHET) 
  call boundary_loss(p)
  CALL localisation(p)

END SUBROUTINE load_part_file
!---------------------------------------------------------------------------
! DESCRIPTION:
!> @brief Checks for each particle if the z position is outside of the local/global simulation space.
!> Depending on the boundary conditions, the leaving particles are sent to the correct neighbouring MPI process
!> or deleted.
!
!> @param[in]  p particles structure
!
!> @author Guillaume Le Bars EPFL/SPC
!---------------------------------------------------------------------------
SUBROUTINE bound(p)

  USE basic, ONLY: zgrid, nz, Zbounds, mpirank, step, leftproc, rightproc, partperiodic
  IMPLICIT NONE
  type(particles), INTENT(INOUT):: p
  INTEGER :: i, rsendnbparts, lsendnbparts, nblostparts
  INTEGER :: receivednbparts, partdiff
  INTEGER, DIMENSION(p%Nploc) :: sendhole
  INTEGER, DIMENSION(p%Nploc) :: losthole
  LOGICAL:: leftcomm, rightcomm
  INTEGER, ALLOCATABLE:: partstoremove(:)

  receivednbparts=0
  nblostparts=0
  rsendnbparts=0
  lsendnbparts=0

  IF (p%Nploc .gt. 0) THEN
  losthole=0
  sendhole=0

  ! We communicate with the left processus
  leftcomm  = leftproc .ne. -1
  ! We communicate with the right processus
  rightcomm = rightproc .ne. -1

  ! Boundary condition at z direction
  !$OMP PARALLEL DO DEFAULT(SHARED)
    DO i=1,p%Nploc
      ! If the particle is to the right of the local simulation space, it is sent to the right MPI process
      IF (p%Z(i) .ge. zgrid(Zbounds(mpirank+1))) THEN
        IF(partperiodic) THEN
          DO WHILE (p%Z(i) .GT. zgrid(nz))
            p%Z(i) = p%Z(i) - zgrid(nz) + zgrid(0)
          END DO
        END IF
        !$OMP CRITICAL (nbparts)
        IF(rightcomm) THEN
            rsendnbparts=rsendnbparts+1
            sendhole(lsendnbparts+rsendnbparts)=i
        ELSE
            nblostparts=nblostparts+1
            losthole(nblostparts)=i
            p%nblost(2)=p%nblost(2)+1
        END IF
        !$OMP END CRITICAL (nbparts)
      ! If the particle is to the left of the local simulation space, it is sent to the left MPI process
      ELSE IF (p%Z(i) .lt. zgrid(Zbounds(mpirank))) THEN
        IF(partperiodic) THEN
          DO WHILE (p%Z(i) .LT. zgrid(0))
            p%Z(i) = p%Z(i) + zgrid(nz) - zgrid(0)
          END DO
        END IF
        !$OMP CRITICAL (nbparts)
        IF(leftcomm) THEN
          ! We send the particle to the left process
          lsendnbparts=lsendnbparts+1
          sendhole(lsendnbparts+rsendnbparts)=-i
        ELSE
          ! we destroy the particle
          nblostparts=nblostparts+1
          losthole(nblostparts)=i
          p%nblost(1)=p%nblost(1)+1
        END IF
        !$OMP END CRITICAL (nbparts)
      END IF
    END DO
  !$OMP END PARALLEL DO
  END IF

    IF(mpisize .gt. 1) THEN
      ! We send the particles leaving the local simulation space to the closest neighbour
      CALL particlescommunication(p, lsendnbparts, rsendnbparts, sendhole, receivednbparts, (/leftproc,rightproc/))
    END IF

    ! If the boundary conditions are not periodic, we delete the corresponding particles
    IF(nblostparts .gt. 0 .and. step .ne. 0) THEN
      DO i=1,nblostparts
          CALL delete_part(p, losthole(i), .false. )
      END DO
      !WRITE(*,'(i8.2,a,i4.2)') nblostparts, " particles lost in z on process: ", mpirank
    END IF

    ! computes if we received less particles than we sent
    partdiff=max(lsendnbparts+rsendnbparts-receivednbparts,0)
    
    IF(nblostparts + partdiff .gt. 0) THEN
      ALLOCATE(partstoremove(nblostparts+partdiff))
      partstoremove(1:partdiff)=abs(sendhole(receivednbparts+1:receivednbparts+partdiff))
      partstoremove(partdiff+1:partdiff+nblostparts)=abs(losthole(1:nblostparts))
      call LSDRADIXSORT(partstoremove,size(partstoremove))
      ! If we received less particles than we sent, or lost particles we fill the remaining holes with the particles from the end of the parts arrays
      DO i=nblostparts+partdiff,1,-1
        CALL move_part(p, p%Nploc, partstoremove(i))
        p%partindex(p%Nploc)=-1
        p%Nploc = p%Nploc-1
      END DO
    END IF
END subroutine bound

   !---------------------------------------------------------------------------
   ! add below the usage of module iiee



   !---------------------------------------------------------------------------
   !> @author
   !> Guillaume Le Bars EPFL/SPC
   !
   ! DESCRIPTION:
   !> @brief Check if a particle is outside the simulation domain and remove it if needed
   !> @param[in]  p particles structure
   !
   !> Last modified on 11/15/2022
   !> Added call to iiee function (see ion_induced_mod.f90) (S.Guinchard) 
   !
   !---------------------------------------------------------------------------
SUBROUTINE boundary_loss(p)
  USE basic, ONLY: rgrid, nr
  Use geometry, ONLY: r_a, geom_weight, dom_weight
  USE iiee
  type(particles), INTENT(INOUT):: p
  INTEGER :: i,j,isup, nblostparts, iend,nbunch
  INTEGER, DIMENSION(p%Nploc) :: losthole
  INTEGER, DIMENSION(16)::idwall
  INTEGER :: nblost(size(p%nblost,1))

  nblostparts=0
  nblost=0
  nbunch=16

  IF (p%Nploc .le. 0) return
  losthole=0

  !$OMP PARALLEL DEFAULT(SHARED), private(i,iend,j,isup,idwall)
      !$OMP DO reduction(+:nblost)
      DO i=1,p%Nploc,nbunch
        ! Avoid segmentation fault caused by accessing non relevant data
        iend=min(i+nbunch-1,p%Nploc)
        ! calculate the weight do determine if a particle is inside the simulation domain.
        call dom_weight(p%Z(i:iend), p%r(i:iend), p%geomweight(i:iend,0),idwall(1:iend-i+1))
        do j=i,iend
          if(p%geomweight(j,0).le.0 .or. p%R(j) .ge. rgrid(nr) .or. p%R(j) .le. rgrid(0)) then
            ! If the particle is outside of the simulation space in the r direction, or if it is outside of the vacuum region it is deleted.
              !$OMP CRITICAL (lostparts)
                nblostparts=nblostparts+1
                losthole(nblostparts)=j
              !$OMP END CRITICAL (lostparts)
                isup=0
                if(p%R(j) .ge. rgrid(nr) .or. idwall(j-i+1) .gt.0) then
                  isup=1
                end if
                nblost(3+isup+idwall(j-i+1))=nblost(3+isup+idwall(j-i+1))+1
          end if
        end do
      END DO
      !$OMP END DO
  !$OMP END PARALLEL
    IF(nblostparts.gt.0) THEN
      p%nblost=nblost+p%nblost
      !call qsort(losthole,p%Nploc,sizeof(losthole(1)),compare_int)
      call LSDRADIXSORT(losthole(1:nblostparts),nblostparts)
      !Write(*,'(a,60i)') "losthole: ", losthole(nblostparts:nblostparts+1)
      IF(p%iiee_id.gt.0) THEN
             CALL ion_induced(p, losthole, partslist(p%iiee_id), nblostparts) 
             !----------------------------------------------------------
             ! CALL ion_induced(p,losthole,partslist(indpelec))
             ! here we call our routine to create electrons out of 
             ! eliminated ions. 
             ! need to define in this file: indpelec (need not to since)
             ! we have the index p%iiee_id
             !----------------------------------------------------------
      END IF 

      DO i=nblostparts,1,-1
          CALL delete_part(p,losthole(i))
      END DO
      
    END IF
  END SUBROUTINE boundary_loss

  !---------------------------------------------------------------------------
   !> @author
   !> Guillaume Le Bars EPFL/SPC
   !
   ! DESCRIPTION:
   !> @brief Compute the grid cell indices for each particle as well as the distance weight Wr, Wz.
   !> @param[in]  p particles structure
   !---------------------------------------------------------------------------
SUBROUTINE localisation(p)
  USE basic, ONLY: rgrid, nr
  Use geometry, ONLY: r_a, geom_weight, dom_weight
  type(particles), INTENT(INOUT):: p
  INTEGER :: i, j,  iend,nbunch
  nbunch=16

  IF (p%Nploc .le. 0) return
    !$OMP PARALLEL DEFAULT(SHARED), private(i,iend,j)
    !$OMP DO
      DO i=1,p%Nploc,nbunch
        ! Avoid segmentation fault by accessing non relevant data
        iend=min(i+nbunch-1,p%Nploc)
        do j=i,iend
          call p_calc_rzindex(p,j)
        end do
        call geom_weight(p%Z(i:iend), p%r(i:iend), p%geomweight(i:iend,:))
      END DO
    !$OMP END DO
    !$OMP END PARALLEL
END SUBROUTINE localisation

subroutine p_calc_rzindex(p,i)
  use basic, only: rgrid,zgrid,invdz,invdr, nnr, nr, nsubr
  integer::i,j,k
  type(particles)::p
  k=0
  do j=1,nsubr
    IF (p%R(i) .GT. rgrid(k) .AND.  p%R(i) .LT. rgrid(k+nnr(j))) THEN
      p%rindex(i)=floor((p%R(i)-rgrid(k))*invdr(j))+k
      exit
    end if
    k=k+nnr(j)
  end do
  !ELSE IF(p%R(i) .GE. rgrid(nnr(1)) .AND.  p%R(i) .LT. rgrid(nnr(1)+nnr(2))) THEN
  !   p%rindex(i)=floor((p%R(i)-rgrid(nnr(1)))*invdr(2))+nnr(1)
  !ELSE IF(p%R(i) .GE. rgrid(nnr(1)+nnr(2)) .AND.  p%R(i) .LT. rgrid(nr)) THEN
  !   p%rindex(i)=floor((p%R(i)-rgrid(nnr(1)+nnr(2)))*invdr(3))+nnr(1)+nnr(2)
  !End if
  p%zindex(i)=floor((p%Z(i)-zgrid(0))*invdz)

end subroutine p_calc_rzindex

SUBROUTINE comp_mag_p(p)
  USE basic, ONLY: zgrid, rgrid, BZ, BR, nz, invdz
  type(particles), INTENT(INOUT):: p
  INTEGER :: i
  Real(kind=db):: WZ,WR
  INTEGER:: j1,j2,j3,j4

  !$OMP PARALLEL DO SIMD DEFAULT(SHARED) Private(J1,J2,J3,J4,WZ,WR)
  DO i=1,p%Nploc
    WZ=(p%Z(i)-zgrid(p%zindex(i)))*invdz;
    WR=(p%R(i)-rgrid(p%rindex(i)))/(rgrid(p%rindex(i)+1)-rgrid(p%rindex(i)));
    J1=(p%rindex(i))*(nz+1) + p%zindex(i)+1
    J2=(p%rindex(i))*(nz+1) + p%zindex(i)+2
    J3=(p%rindex(i)+1)*(nz+1)+p%zindex(i)+1
    J4=(p%rindex(i)+1)*(nz+1)+p%zindex(i)+2
    
    ! Interpolation for magnetic field
    p%BZ(i)=(1-WZ)*(1-WR)*Bz(J4) &
    &   +WZ*(1-WR)*Bz(J3)    &
    &   +(1-WZ)*WR*Bz(J2)    &
    &   +WZ*WR*Bz(J1)
    p%BR(i)=(1-WZ)*(1-WR)*Br(J4) &
    &   +WZ*(1-WR)*Br(J3)    &
    &   +(1-WZ)*WR*Br(J2)    &
    &   +WZ*WR*Br(J1)
  END DO
  !$OMP END PARALLEL DO SIMD
end subroutine comp_mag_p

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Routine used to compute the lorentz factor \f$\gamma\f$ in the classical simulations.
!> This routine systematically returns 1.0 to treat the system according to classical dynamic.
!
!> @param[out] gamma the lorentz factor \f$\gamma\f$
!> @param[in]  UZ \f$\gamma\beta_z=\gamma v_z/c\f$ the normalized particle longitudinal velocity
!> @param[in]  UR \f$\gamma\beta_r=\gamma v_r/c\f$ the normalized particle radial velocity
!> @param[in]  UTHET \f$\gamma\beta_\theta=\gamma v_\theta/c\f$ the normalized particle azimuthal velocity
!---------------------------------------------------------------------------
REAL(kind=db) FUNCTION gamma_classical(UZ, UR, UTHET)
#if __INTEL_COMPILER > 1700
!$OMP declare simd(gamma_classical)
#endif
      REAL(kind=db), INTENT(IN):: UR,UZ,UTHET
      gamma_classical=1.0
END FUNCTION gamma_classical
!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief Routine used to compute the lorentz factor \f$\gamma\f$ in the relativistic simulations.
!> This routine computes the Lorentz factor \f$\gamma=\sqrt{1+\mathbf{\gamma\beta}^2}\f$
!
!> @param[out] gamma the lorentz factor \f$\gamma\f$
!> @param[in]  UZ \f$\gamma\beta_z=\gamma v_z/c\f$ the normalized particle longitudinal velocity
!> @param[in]  UR \f$\gamma\beta_r=\gamma v_r/c\f$ the normalized particle radial velocity
!> @param[in]  UTHET \f$\gamma\beta_\theta=\gamma v_\theta/c\f$ the normalized particle azimuthal velocity
!---------------------------------------------------------------------------
REAL(kind=db) FUNCTION gamma_relativistic(UZ, UR, UTHET)
#if __INTEL_COMPILER > 1700
!$OMP declare simd(gamma_relativistic)
#endif
      REAL(kind=db), INTENT(IN):: UR,UZ,UTHET
      gamma_relativistic=sqrt(1+UZ**2+UR**2+UTHET**2)
END FUNCTION gamma_relativistic

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief General routine to compute the velocities at time t+1.
!> This routine allows to treat the classical and relativistic case efficiently from a numerical standpoint,
!> by using a pointer to the routine computing gamma. This avoid the nlclassical flag check on each particle.
!
!> @param[in]  p The particles structure being updated
!---------------------------------------------------------------------------
SUBROUTINE comp_velocity(p)
!
!   Computes the new velocity of the particles due to Lorentz force
!
  USE basic, ONLY : nlclassical
  type(particles), INTENT(INOUT):: p
  ! Store old Velocities
  CALL swappointer(p%UZold,    p%UZ)
  CALL swappointer(p%URold,    p%UR)
  CALL swappointer(p%UTHETold, p%UTHET)
  CALL swappointer(p%Gammaold, p%Gamma)

  IF (nlclassical) THEN
    CALL comp_velocity_fun(p, gamma_classical)
  ELSE
    CALL comp_velocity_fun(p, gamma_relativistic)
  END IF

END SUBROUTINE comp_velocity

!---------------------------------------------------------------------------
!> @author
!> Patryk Kaminski EPFL/SPC
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief Routine called by comp_velocity to compute the velocities at time t+1.
!> This routine allows to treat the classical and relativistic case efficiently from a numerical standpoint,
!> by using the routine computing gamma as an input. This avoid the nlclassical flag check on each particle.
!
!> @param[in] gamma the function used to compute the value of the lorentz factor \f$\gamma\f$
!> @param[in]  p The particles structure being updated
!---------------------------------------------------------------------------
SUBROUTINE comp_velocity_fun(p, gammafun)
!
!   Computes the new velocity of the particles due to Lorentz force
!
  USE basic, ONLY : bnorm, dt, tnorm
  procedure(gamma)::gammafun
  type(particles), INTENT(INOUT):: p
  REAL(kind=db) :: tau
  REAL(kind=db):: BRZ, BRR, ZBR, ZBZ, ZPR, ZPZ, ZPTHET, SQR, ZBZ2, ZBR2
  INTEGER:: J1, J2, J3, J4
  INTEGER:: i
 ! Normalized time increment
  !tau=p%qmRatio*bnorm*tnorm*0.5*dt/tnorm
  tau=p%qmRatio*bnorm*0.5*dt*tnorm
  IF (p%Nploc .NE. 0) THEN
    !$OMP PARALLEL DO SIMD DEFAULT(SHARED) PRIVATE(J1,J2,J3,J4,BRZ, BRR, ZBR, ZBZ, ZPR, ZPZ, ZPTHET, SQR, ZBZ2, ZBR2)
    DO i=1,p%Nploc
    ! First half of electric pulse
       p%UZ(i)=p%UZold(i)+p%Ez(i)*tau
       p%UR(i)=p%URold(i)+p%ER(i)*tau
       p%Gamma(i)=gammafun(p%UZ(i), p%UR(i), p%UTHETold(i))

    ! Rotation along magnetic field
       ZBZ=tau*p%BZ(i)/p%Gamma(i)
       ZBR=tau*p%BR(i)/p%Gamma(i)
       ZPZ=p%UZ(i)-ZBR*p%UTHETold(i)                     !u'_{z}
       ZPR=p%UR(i)+ZBZ*p%UTHETold(i)                     !u'_{r}
       ZPTHET=p%UTHETold(i)+(ZBR*p%UZ(i)-ZBZ*p%UR(i))          !u'_{theta}
       SQR=1+ZBZ*ZBZ+ZBR*ZBR
       ZBZ2=2*ZBZ/SQR
       ZBR2=2*ZBR/SQR
       p%UZ(i)=p%UZ(i)-ZBR2*ZPTHET                    !u+_{z}
       p%UR(i)=p%UR(i)+ZBZ2*ZPTHET                   !u+_{r}
       p%UTHET(i)=p%UTHETold(i)+(ZBR2*ZPZ-ZBZ2*ZPR)      !u+_{theta}

    ! Second half of acceleration
       p%UZ(i)=p%UZ(i)+p%EZ(i)*tau
       p%UR(i)=p%UR(i)+p%ER(i)*tau
       !p%ur(i)=0.001
    ! Final computation of the Lorentz factor
       p%Gamma(i)=gammafun(p%UZ(i), p%UR(i), p%UTHET(i))
    END DO
    !$OMP END PARALLEL DO SIMD
  END IF
  p%collected=.false.
END SUBROUTINE comp_velocity_fun

!---------------------------------------------------------------------------
!> @author
!> Patryk Kaminski EPFL/SPC
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Computes the particles position at time t+1
!> This routine computes the particles position at time t+1 according to the Bunemann algorithm.
!
!> @param[in]  p The particles structure being updated
!---------------------------------------------------------------------------
SUBROUTINE push(p)
    Use basic, ONLY: dt, tnorm
    type(particles), INTENT(INOUT):: p
    REAL(kind=db):: XP, YP, COSA, SINA, U1, U2, ALPHA
    INTEGER :: i

    IF (p%Nploc .NE. 0) THEN
    !$OMP PARALLEL DO SIMD DEFAULT(SHARED) PRIVATE(XP, YP, COSA, SINA, U1, U2, ALPHA)
    DO i=1,p%Nploc
! Local Cartesian coordinates
          XP=p%R(i)+dt*p%UR(i)/p%Gamma(i)
          YP=dt*p%UTHET(i)/p%Gamma(i)

        ! Conversion to cylindrical coordiantes
          p%Z(i)=p%Z(i)+dt*p%UZ(i)/p%Gamma(i)
          p%R(i)=sqrt(XP**2+YP**2)

        ! Computation of the rotation angle
          IF (p%R(i) .EQ. 0) THEN
            COSA=1
            SINA=0
            ALPHA=0
          ELSE
            COSA=XP/p%R(i)
            SINA=YP/p%R(i)
            ALPHA=asin(SINA)
          END IF
        ! New azimuthal position
          p%THET(i)=MOD(p%THET(i)+ALPHA,2*pi)

        ! Velocity in rotated reference frame
          U1=COSA*p%UR(i)+SINA*p%UTHET(i)
          U2=-SINA*p%UR(i)+COSA*p%UTHET(i)

          p%UR(i)=U1
          p%UTHET(i)=U2

        END DO
        !$OMP END PARALLEL DO SIMD
    END IF
    p%collected=.false.
END SUBROUTINE push
!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Computes several diagnostic quantities
!> This routine computes the total kinetic and electric potential energy.
!> It keeps track of the reference energy and the number of particle per mpi node. 
!
!---------------------------------------------------------------------------
SUBROUTINE partdiagnostics
!
!  Compute energies
!
    USE constants, ONLY: vlight
    USE basic, ONLY: phinorm, cstep, nlclassical, ierr, step, nlend,&
    & itparts, nbspecies

    INTEGER:: i,j

  ! Reset the quantities
    ekin=0
    epot=0
    etot=0


  ! Computation of the kinetic and potential energy as well as  fluid velocities and density
    !$OMP PARALLEL DO REDUCTION(+:epot, ekin) DEFAULT(SHARED), PRIVATE(i,j) 
    Do j=1,nbspecies
      if(.not. partslist(j)%is_field) CYCLE
      DO i=1,partslist(j)%Nploc

!        Potential energy
         epot=epot+(partslist(j)%pot(i)+partslist(j)%potxt(i))*partslist(j)%q*partslist(j)%weight
         ! Kinetic energy
         IF(.not. nlclassical) THEN
          ekin=ekin+(0.5*(partslist(j)%Gammaold(i)+partslist(j)%Gamma(i))-1)*partslist(j)%m*partslist(j)%weight
         ELSE
          ekin=ekin+0.5*( partslist(j)%UR(i)*partslist(j)%URold(i)    &
                      & + partslist(j)%UZ(i)*partslist(j)%UZold(i)    &
                      & + partslist(j)%UTHET(i)*partslist(j)%UTHETold(i) )*partslist(j)%m*partslist(j)%weight
         END IF
      END DO
    END DO
    !$OMP END PARALLEL DO
    epot=epot*phinorm*0.5
    ekin=ekin*vlight**2

    !  Shift to Etot at cstep=1 (not valable yet at cstep=0!)
    IF(cstep.EQ. 1) THEN
    ! Compute the local total energy
       loc_etot0 = epot+ekin
       etot0=0
    END IF
    !etot=loc_etot0
    ! Compute the total energy
    etot=epot+ekin
    Energies=(/ekin,epot,etot,loc_etot0/)
  ! The computed energy is sent to the root process
    IF(mpisize .gt.1) THEN
      IF(mpirank .eq.0 ) THEN
          CALL MPI_REDUCE(MPI_IN_PLACE, Energies, 4, db_type, db_sum_op, &
          & 0, MPI_COMM_WORLD, ierr)
          etot0=etot0+Energies(4)
          ekin=Energies(1)
          epot=Energies(2)
          etot=Energies(3)
      ELSE
          CALL MPI_REDUCE(Energies, Energies, 4, db_type, db_sum_op, &
          & 0, MPI_COMM_WORLD, ierr)
      END IF
    ELSE
      etot0=etot0+loc_etot0
    END IF
    loc_etot0=0

  ! Send the local number of particles on each node to the root process
    IF(mpisize .gt. 1) THEN
     Nplocs_all(mpirank)=partslist(1)%Nploc
     IF(mpirank .eq.0 ) THEN
        CALL MPI_gather(MPI_IN_PLACE, 1, MPI_INTEGER, Nplocs_all, 1, MPI_INTEGER,&
        & 0, MPI_COMM_WORLD, ierr)
        !CALL MPI_REDUCE(MPI_IN_PLACE,partslist(1)%nudcol,3,db_type,db_sum_op,0,MPI_COMM_WORLD,ierr)
        partslist(1)%Nptot=sum(Nplocs_all)
       !partslist(1)%nudcol=partslist(1)%nudcol/partslist(1)%Nptot
      ELSE
        CALL MPI_gather(Nplocs_all(mpirank), 1, MPI_INTEGER, Nplocs_all, 1, MPI_INTEGER,&
        & 0, MPI_COMM_WORLD, ierr)
        !CALL MPI_REDUCE(partslist(1)%nudcol,partslist(1)%nudcol,3,db_type,db_sum_op,0,MPI_COMM_WORLD,ierr)
      END IF
    ELSE
      partslist(1)%Nptot=partslist(1)%Nploc
    END IF
  end subroutine partdiagnostics

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief Collect the particles positions and velocities on the root process.
!> If the collection has already been performed at this time step, the routine does nothing.
!
!---------------------------------------------------------------------------
  SUBROUTINE collectparts(p)
    USE basic, ONLY: mpirank, mpisize, ierr
    type(particles), INTENT(INOUT):: p
    INTEGER, DIMENSION(:), ALLOCATABLE :: displs, Nploc
    INTEGER:: i
    INTEGER:: particles_type(mpisize-1)     !< Stores the MPI data type used for particles gathering on node 0 and broadcast from node 0
    INTEGER :: part_requests(mpisize-1)
    INTEGER:: stats(MPI_STATUS_SIZE,mpisize-1)

    part_requests=MPI_REQUEST_NULL

    particles_type=MPI_DATATYPE_NULL
    IF(p%collected) RETURN ! exit subroutine if particles have already been collected during this time step

    ALLOCATE(Nploc(0:mpisize-1))
    ALLOCATE(displs(0:mpisize-1))
    displs=0

    Nploc(mpirank)=p%Nploc
    CALL MPI_Allgather(MPI_IN_PLACE, 1, MPI_INTEGER, Nploc, 1, MPI_INTEGER,&
    & MPI_COMM_WORLD, ierr)
    
    p%Nptot=sum(Nploc)
    IF(p%Nptot .eq. 0 ) THEN
      p%partindex=-1
      p%collected=.true.
      RETURN
    END IF

    Do i=1,mpisize-1
      displs(i)=displs(i-1)+Nploc(i-1)
    END DO
    IF(mpirank.eq.0 .and. p%Nptot .gt. size(p%R,1)) THEN
     CALL change_parts_allocation(p,max(p%Nptot-size(P%R,1),floor(0.5*size(P%R,1))))
   END IF

    

    IF(mpirank .ne. 0) THEN
      if(Nploc(mpirank) .gt. 0) THEN
        Call init_particles_gather_mpi(p,1,Nploc(mpirank),particles_type(mpirank))
        ! Send Particles informations to root process
        CALL MPI_SEND(p, 1, particles_type(mpirank), 0, partsgather_tag, MPI_COMM_WORLD, ierr)
        CALL MPI_TYPE_FREE(particles_type(mpirank),ierr)
      END IF
    ELSE
  ! Receive particle information from all processes
      DO i=1,mpisize-1
        if(Nploc(i) .lt. 1) cycle
        Call init_particles_gather_mpi(p,displs(i)+1,Nploc(i),particles_type(i))
        CALL MPI_IRECV(p,1,particles_type(i),i,partsgather_tag,MPI_COMM_WORLD, part_requests(i), ierr)
      END DO
      CALL MPI_WAITALL(mpisize-1,part_requests, stats, ierr)
      p%partindex(sum(Nploc)+1:)=-1
      Do i=1,mpisize-1
        if(Nploc(i) .lt. 1) cycle
        CALL MPI_TYPE_FREE(particles_type(i),ierr)
      END DO
    END IF
    p%collected=.TRUE.
  END SUBROUTINE collectparts
!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief Computes the velocities at time t-1/2 delta t to keep the second order precision in time on the velocity.
!> This should only be used at particle initialisation time, ot in the case of a restart.
!
!---------------------------------------------------------------------------
  SUBROUTINE adapt_vinit(p)
    !!   Computes the velocity at time -dt/2 from velocities computed at time 0
    !
      USE basic, ONLY : bnorm, dt, tnorm, nlclassical, phinorm, distribtype, vnorm
      type(particles), INTENT(INOUT):: p
      REAL(kind=db) :: tau, BRZ, BRR, ZBR, ZBZ, ZPR, ZPZ, ZPTHET, &
      &           SQR, Vperp, v2
      INTEGER :: J1, J2, J3, J4, i
      REAL(kind=db), DIMENSION(:), ALLOCATABLE :: VZ, VR, VTHET

    ! In case Davidson distribution is used the longitudinal and radial velocities are adapted to take into account the
    ! electric potential.
      IF(distribtype .EQ. 2 .OR. distribtype .EQ. 3 .OR. distribtype .EQ. 4 .or. p%Davidson) THEN
        ALLOCATE(VR(p%Nploc),VZ(p%Nploc),VTHET(p%Nploc))
        CALL loduni(7,VZ)
        VZ=VZ*2*pi
        VTHET=p%UTHET/p%Gamma*vnorm
        DO i=1,p%Nploc
          Vperp=sqrt(MAX(2*p%H0/p%m-2*p%qmRatio*p%pot(i)*phinorm-VTHET(i)**2,0.0_db))
          VR(i)=Vperp*sin(VZ(i))
          VZ(i)=Vperp*cos(VZ(i))
          IF(nlclassical) THEN
            p%Gamma(i)=1
          ELSE
            v2=VR(i)**2+VZ(i)**2+VTHET(i)**2
            p%Gamma(i)=sqrt(1/(1-v2/vnorm**2))
          END IF
          p%UR(i)=p%Gamma(i)*VR(i)/vnorm
          p%UZ(i)=p%Gamma(i)*VZ(i)/vnorm
          p%UTHET(i)=p%Gamma(i)*VTHET(i)/vnorm
        END DO
        DEALLOCATE(VR,VZ,VTHET)
      END IF

      ! Normalised time increment
      !tau=-omegac/2/omegap*dt/tnorm
      tau=-p%qmRatio*bnorm*0.5*dt*tnorm
      ! Store old Velocities
      CALL swappointer(p%UZold, p%UZ)
      CALL swappointer(p%URold, p%UR)
      CALL swappointer(p%UTHETold, p%UTHET)
      CALL swappointer(p%Gammaold, p%Gamma)

      IF (p%Nploc .NE. 0) THEN
      !$OMP PARALLEL DO SIMD DEFAULT(SHARED) PRIVATE(J1,J2,J3,J4,BRZ, BRR, ZBR, ZBZ, ZPR, ZPZ, ZPTHET, SQR)
        DO i=1,p%Nploc

        ! Half inverse Rotation along magnetic field
           ZBZ=tau*p%BZ(i)/p%Gammaold(i)
           ZBR=tau*p%BR(i)/p%Gammaold(i)
           SQR=1+ZBZ*ZBZ+ZBR*ZBR
           ZPZ=(p%UZold(i)-ZBR*p%UTHETold(i))/SQR                     !u-_{z}
           ZPR=(p%URold(i)+ZBZ*p%UTHETold(i))/SQR                     !u-_{r}
           ZPTHET=p%UTHETold(i)+(ZBR*p%UZold(i)-ZBZ*p%URold(i))/SQR          !u-_{theta}
           p%UZ(i)=ZPZ
           p%UR(i)=ZPR
           p%UTHET(i)=ZPTHET

        ! half of decceleration
           p%UZ(i)=p%UZ(i)+p%Ez(i)*tau
           p%UR(i)=p%UR(i)+p%Er(i)*tau

           IF(.not. nlclassical) THEN
            p%Gamma(i)=sqrt(1+p%UZ(i)**2+p%UR(i)**2+p%UTHET(i)**2)
           END IF

        END DO
        !$OMP END PARALLEL DO SIMD
      END IF

    END SUBROUTINE adapt_vinit

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief Calculates the number of paticles per column of the spatial grid ( at fixed axial cell position)
!> This facilitate the computation of the axial grid limits for each MPI worker
!
!---------------------------------------------------------------------------
SUBROUTINE calcnbperz(p,nbperz)
  USE basic, only: nz
  IMPLICIT NONE
  type(particles):: p
  INTEGER, INTENT(INOUT):: nbperz(0:)
  Integer::i, zindex 
  
  nbperz=0
  !$OMP PARALLEL DO DEFAULT(SHARED) reduction(+:nbperz), private(zindex,i)
  Do i=1,p%Nploc
    ! we make sure zindex is in [0, nz-1] to avoid segmentation faults
    zindex=min(nz-1,max(p%zindex(i),0)) 
    nbperz(zindex)=nbperz(zindex)+1
  END DO
  !$OMP END PARALLEL DO
END SUBROUTINE calcnbperz

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief In the case of MPI parallelism, computes the axial limits of the local domain.
!---------------------------------------------------------------------------
  SUBROUTINE calc_Zbounds(p, Zbounds, norder)
  ! Computes the start and end indices for the Z boundaries on local processus
  ! Computes the particle indices from initial particle loading vector, that stay in current process
    USE basic, ONLY: nz, cstep, mpirank, mpisize,step
    USE mpihelper
    TYPE(particles), INTENT(INOUT):: p
    INTEGER:: Zbounds(0:)
    INTEGER:: norder(2)
    INTEGER:: old_Zbounds(0:size(Zbounds,1)-1)
    INTEGER:: k, i, nbparts
    REAL(kind=db):: idealnbpartsperproc
    INTEGER, DIMENSION(0:nz-1):: partspercol ! Vector containing the number of particles between zgrid(n) and zgrid(n+1)
    INTEGER:: Zmin, Zmax ! Minimum and maximum indices of particles in Z direction
    INTEGER:: Zperproc, ierr, remparts
    CHARACTER(12)::fmt

    ! calculatese the axial disstibution integrated along the radial direction
    call calcnbperz(p,partspercol)

    ! gather this data on all nodes
    if(step .gt. 0 .and. mpisize .gt. 1) THEN
      old_Zbounds=Zbounds
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, partspercol, nz, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    END IF

    ! estimate the ideal number of particle per MPI worker
    idealnbpartsperproc = p%Nptot/mpisize

    ! find the start and end indices where particles are present
    Zmin=0
    Zmax=nz-1
    Do k=0,nz-1
      if(partspercol(k) .gt.0) then
        Zmin=k
        exit
      end if
    end do
    Do k=nz-1,0,-1
      if(partspercol(k) .gt.0) then
        Zmax=k
        exit
      end if
    end do 
    !Zmin=findloc(partspercol.gt.0,.true.,1,kind=kind(Zmin))-1
    !Zmax=findloc(partspercol.gt.0,.true.,1,kind=kind(Zmax),BACK=.true.)-1

    ! Find naive axial limits assuming uniform axial distribution
    IF(Zmax .le. 0) Zmax=nz-1
    IF(Zmin .gt. nz) Zmin=0
    Zperproc=(Zmax-Zmin)/mpisize

   
    IF (Zperproc .lt. 1 .or. cstep .eq. 0) THEN
      !! No particles are present initially
      Zperproc=nz/mpisize
      Zmin=0
      ! Define boundaries using naive guess on start or restart (allow to start with 0 parts)
      DO k=1,mpisize-1
        IF(k .lt. mpisize-1-MODULO(Zmax-Zmin,mpisize)) THEN
          Zbounds(k)=Zmin+k*Zperproc-1
        ELSE
          Zbounds(k)=Zmin+k*Zperproc-1+k-mpisize+2+MODULO(Zmax-Zmin,mpisize)
        END IF
      END DO

    ELSE
      i=0
      ! Define axial boundaries using the axial distribution information.
      ! the subdomains are not equal
      remparts=p%Nptot
      DO k=1,mpisize-1
        nbparts=0
        DO WHILE(nbparts<0.98*idealnbpartsperproc .and. i .lt. Zmax .and. (nbparts+partspercol(i)).lt.1.25*idealnbpartsperproc)
          nbparts=nbparts+partspercol(i)
          i=i+1 
        END DO
        remparts=remparts-nbparts
        Zbounds(k)=i
      END DO
    END IF

    IF(step .gt. 0 .and. mpirank .eq. 0) THEN
      Do i=1,mpisize-1
        !We check that the new limits will not exceed the old limits of the left and right process
        ! This avoids particle communications with process >mpirank+2 and < mpirank-2
        ! However this should converge over time
        IF(Zbounds(i) .lt. old_Zbounds(i-1)) Zbounds(i)=old_Zbounds(i-1)
        if(Zbounds(i) .gt. old_Zbounds(i+1))Zbounds(i)=old_Zbounds(i+1)
        ! If a process would have an axial domain shoter than axial norder, we revert to the old boundaries.
        IF((Zbounds(i)-Zbounds(i-1)).lt. norder(1) .or. (Zbounds(i+1)-Zbounds(i)).lt. norder(1)) THEN
          Zbounds=old_Zbounds
          EXIT
        END IF
      END DO
    END IF

    ! send the new boundaries to all the workers
    CALL MPI_Bcast(Zbounds,mpisize+1,MPI_INTEGER,0,MPI_COMM_WORLD, ierr)
    DO k=0,mpisize-1
      Nplocs_all(k)=SUM(partspercol(Zbounds(k):Zbounds(k+1)-1))
    END DO

    if(mpirank .eq. 0) THEN
      WRITE(fmt,'(a,i3,a)')"(a,",mpisize+1, "i5)"
      WRITE(*,fmt) "Zbounds: ", Zbounds
      WRITE(fmt,'(a,i3,a)')"(a,",mpisize, "i8)"
      WRITE(*,fmt) "Nplocs: ", Nplocs_all
    END IF
    
  END SUBROUTINE calc_Zbounds

  !---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!> @brief After a restart keep only the particles in the local domain of the current MPI worker
!---------------------------------------------------------------------------
  SUBROUTINE keep_mpi_self_parts(p,Zbounds)
    TYPE(particles),INTENT(INOUT):: p
    INTEGER,INTENT(in)::Zbounds(0:)
    INTEGER :: i, partstart, old_sum,ierr
    partstart=1
    p%Nploc=0
    Do i=1,p%Nptot
      IF(p%Zindex(i).ge.Zbounds(mpirank).and.p%Zindex(i).lt.Zbounds(mpirank+1)) THEN
        p%Nploc=p%Nploc+1
        CALL move_part(p,i,p%Nploc)
      END IF
    END DO
    old_sum=p%Nptot
    CALL MPI_REDUCE(p%Nploc, p%Nptot,1,MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    IF(p%Nptot .ne. old_sum) THEN
       WRITE(*,*) "Error in particle distribution kept: ", p%Nptot, "/",old_sum
       !call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
       !stop
    END IF
  END SUBROUTINE keep_mpi_self_parts

!_______________________________________________________________________________
!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Manage the particle communication between neighbours.
!> This routine is responsible to receive the incoming particles from the MPI neighbours and to send its outgoing
!> particles to these neighbours
!
!> @param [in] lsendnbparts number of particles to send to the left neighbour (mpirank-1)
!> @param [in] rsendnbparts number of particles to send to the right neighbour (mpirank+1)
!> @param [in] sendholes array containing the indices of the particle leaving the local domain in ascending order. If the index is positive, the particle goes to the right neigbour, and to the left neighbour if the index is negative
!---------------------------------------------------------------------------
  SUBROUTINE particlescommunication(p, lsendnbparts, rsendnbparts, sendholes, receivednbparts, procs)
    USE mpihelper, ONLY: particle_type
#ifdef _DEBUG
    USE basic, ONLY: step
#endif
    type(particles), INTENT(INOUT):: p
    INTEGER, INTENT(in)    :: lsendnbparts, rsendnbparts
    INTEGER, INTENT(out)   :: receivednbparts
    INTEGER, INTENT(in) :: sendholes(:)
    INTEGER, INTENT(in) :: procs(2)
    INTEGER,  ASYNCHRONOUS :: rrecvnbparts=0, lrecvnbparts=0
    INTEGER, ASYNCHRONOUS :: sendrequest(2), recvrequest(2)
    INTEGER, ASYNCHRONOUS :: sendstatus(MPI_STATUS_SIZE,2), recvstatus(MPI_STATUS_SIZE,2)
    TYPE(particle), ALLOCATABLE :: rrecvpartbuff(:), lrecvpartbuff(:), rsendpartbuff(:), lsendpartbuff(:) ! buffers to send and receive particle from left and right processes

    INTEGER :: lsentnbparts, rsentnbparts
    INTEGER :: lreceivednbparts, rreceivednbparts, ierr

    lsentnbparts=lsendnbparts
    rsentnbparts=rsendnbparts

    sendrequest=MPI_REQUEST_NULL
    recvrequest=MPI_REQUEST_NULL
    lrecvnbparts=0
    rrecvnbparts=0

  ! Send and receive the number of particles to exchange
    CALL MPI_IRECV(lrecvnbparts, 1, MPI_INTEGER, procs(1),  nbpartsexchange_tag, MPI_COMM_WORLD, recvrequest(1), ierr)
    CALL MPI_IRECV(rrecvnbparts, 1, MPI_INTEGER, procs(2), nbpartsexchange_tag, MPI_COMM_WORLD, recvrequest(2), ierr)
    CALL MPI_ISEND(lsentnbparts, 1, MPI_INTEGER, procs(1),  nbpartsexchange_tag, MPI_COMM_WORLD, sendrequest(1), ierr)
    CALL MPI_ISEND(rsentnbparts, 1, MPI_INTEGER, procs(2), nbpartsexchange_tag, MPI_COMM_WORLD, sendrequest(2), ierr)

    CALL MPI_Waitall(2,recvrequest(1:2), recvstatus(:,1:2), ierr)

    recvrequest=MPI_REQUEST_NULL

    lreceivednbparts=lrecvnbparts
    rreceivednbparts=rrecvnbparts

  ! Re/allocate enough memory to store the incoming particles
    ALLOCATE(rrecvpartbuff(rreceivednbparts))
    ALLOCATE(lrecvpartbuff(lreceivednbparts))

  ! Receive particles from left and right processes to the corresponding buffers
    IF ( lrecvnbparts .gt. 0) THEN
      CALL MPI_IRECV(lrecvpartbuff, lreceivednbparts, particle_type, procs(1),  partsexchange_tag, MPI_COMM_WORLD, recvrequest(1), ierr)
    END IF
    IF( rrecvnbparts .gt. 0) THEN
      CALL MPI_IRECV(rrecvpartbuff, rreceivednbparts, particle_type, procs(2), partsexchange_tag, MPI_COMM_WORLD, recvrequest(2), ierr)
    END IF

    ALLOCATE(rsendpartbuff(rsendnbparts))
    ALLOCATE(lsendpartbuff(lsendnbparts))

  ! Copy the leaving particles to the corresponding send buffers
    IF ( (lsendnbparts + rsendnbparts) .gt. 0) THEN
      CALL AddPartSendBuffers(p, lsendnbparts, rsendnbparts, sendholes, lsendpartbuff, rsendpartbuff)
    END IF

    CALL MPI_Waitall(2,sendrequest(1:2), sendstatus(:,1:2), ierr)

  ! Send the particles to the left and right neighbours
    IF( lsendnbparts .gt. 0) THEN
        CALL MPI_ISEND(lsendpartbuff, lsendnbparts, particle_type, procs(1),  partsexchange_tag, MPI_COMM_WORLD, sendrequest(1), ierr)
#ifdef _DEBUG
        !WRITE(*,*)"snding ", lsendnbparts , " to  left  at step: ",step
#endif
    END IF
    IF( rsendnbparts .gt. 0) THEN
        CALL MPI_ISEND(rsendpartbuff, rsendnbparts, particle_type, procs(2), partsexchange_tag, MPI_COMM_WORLD, sendrequest(2), ierr)
#ifdef _DEBUG
        !WRITE(*,*)"snding ", rsendnbparts , " to  right at step: ",step
#endif
    END IF

  ! Receive the incoming parts in the receive buffers
    IF ( lreceivednbparts .gt. 0) THEN
      CALL MPI_Wait(recvrequest(1), recvstatus(:,1), ierr)
      IF(ierr .ne. MPI_SUCCESS) THEN
        WRITE(*,*) "Error in particle reception on proc:", mpirank, " error code:", ierr, "status:", recvstatus(:,1)
        CALL MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      END IF
#ifdef _DEBUG
      !WRITE(*,*)"rcvd ", lreceivednbparts , " from left  at step: ",step
#endif
    END IF
    IF ( rreceivednbparts .gt. 0) THEN
      CALL MPI_Wait(recvrequest(2), recvstatus(:,2), ierr)
      IF(ierr .ne. MPI_SUCCESS) THEN
        WRITE(*,*) "Error in particle reception on proc:", mpirank, " error code:", ierr, "status:", recvstatus(:,2)
        CALL MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      END IF
#ifdef _DEBUG
      !WRITE(*,*)"rcvd ", rreceivednbparts , " from right at step: ",step
#endif
    END IF

    receivednbparts=rreceivednbparts+lreceivednbparts

    IF(p%Nploc+receivednbparts-lsendnbparts-rsendnbparts .gt. size(p%R,1)) THEN
      CALL change_parts_allocation(p,receivednbparts)
    END IF

  ! Copy the incoming particles from the receive buffers to the simulation parts variable
    CALL Addincomingparts(p, rreceivednbparts, lreceivednbparts, lsendnbparts+rsendnbparts, &
    & sendholes, lrecvpartbuff, rrecvpartbuff)

  ! Wait for the outgoing particles to be fully received by the neighbours
    IF( lsendnbparts .gt. 0) THEN
      CALL MPI_Wait(sendrequest(1), sendstatus(:,1), ierr)
#ifdef _DEBUG
      !WRITE(*,*)"sent ", lsentnbparts , " to  left  at step: ",step
#endif
    END IF
    IF( rsendnbparts .gt. 0) THEN
      CALL MPI_Wait(sendrequest(2), sendstatus(:,2), ierr)
#ifdef _DEBUG
      !WRITE(*,*)"sent ", rsentnbparts , " to  right at step: ",step
#endif
    END IF
!
!
  END SUBROUTINE particlescommunication

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Copy the particles from the receive buffers to the local simulation variable parts.
!> The incoming particles will first be stored in the holes left by the outgoing particles, then they
!> will be added at the end of the parts variable
!
!> @param [in] rrecvnbparts number of particles received from the right neighbour (mpirank+1)
!> @param [in] lrecvnbparts number of particles received from the left neighbour (mpirank-1)
!> @param [in] sendnbparts total number of particles having left the local domain
!> @param [in] sendholes array containing the indices of the particle having left the local domain in ascending order.
!---------------------------------------------------------------------------
  SUBROUTINE Addincomingparts(p, rrecvnbparts, lrecvnbparts, sendnbparts, sendholes,lrecvpartbuff, rrecvpartbuff)
!
    USE mpihelper
    TYPE(particles), INTENT(INOUT):: p
    INTEGER, INTENT(in)    :: rrecvnbparts, lrecvnbparts, sendnbparts
    INTEGER, INTENT(in) :: sendholes(:)
    TYPE(particle), INTENT(IN) :: rrecvpartbuff(:), lrecvpartbuff(:) 
    INTEGER k,partpos

  

  ! First import the particles coming from the right
    IF(rrecvnbparts .gt. 0) THEN
      Do k=1,rrecvnbparts
        IF(k .le. sendnbparts) THEN
        ! Fill the holes left by sent parts
          partpos=abs(sendholes(k))
        ELSE
        ! Add at the end of parts and keep track of number of parts
          p%Nploc=p%Nploc+1
          partpos=p%Nploc
        END IF
        CALL Insertincomingpart(p, rrecvpartbuff(k), partpos)
      END DO
    END IF

  ! Then import the particles coming from the left
    IF(lrecvnbparts .gt. 0) THEN
      Do k=1,lrecvnbparts
        IF(k+rrecvnbparts .le. sendnbparts) THEN
          ! Fill the holes left by sent parts
          partpos=abs(sendholes(k+rrecvnbparts))
        ELSE
          ! Add at the end of parts and keep track of number of parts
          p%Nploc=p%Nploc+1
          partpos=p%Nploc
        END IF
        CALL Insertincomingpart(p, lrecvpartbuff(k), partpos)
      END DO
    END IF
!
  END SUBROUTINE Addincomingparts


!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Copy the particles from the local parts variable to the left and right send buffers.
!
!> @param [in] lsendnbparts number of particles to send to the left neighbour (mpirank-1)
!> @param [in] rsendnbparts number of particles to send to the right neighbour (mpirank+1)
!> @param [in] sendholes array containing the indices of the particle leaving the local domain in ascending order. If the index is positive, the particle goes to the right neigbour, and to the left neighbour if the index is negative
!---------------------------------------------------------------------------
  SUBROUTINE AddPartSendBuffers(p, lsendnbparts, rsendnbparts, sendholes, lsendpartbuff, rsendpartbuff)
!
    USE mpihelper
    TYPE(particles), INTENT(INOUT):: p
    INTEGER, INTENT(in)    :: lsendnbparts, rsendnbparts
    INTEGER, INTENT(in) :: sendholes(:)
    TYPE(particle), INTENT(OUT) :: rsendpartbuff(:), lsendpartbuff(:) 
    INTEGER:: partpos, k
    INTEGER:: lsendpos, rsendpos
    lsendpos=0
    rsendpos=0

  ! Loop over the outgoing particles and fill the correct send buffer
    Do k=lsendnbparts+rsendnbparts,1,-1
      partpos=abs(sendholes(k))
      IF(sendholes(k) .GT. 0) THEN
        rsendpos=rsendpos+1
        CALL Insertsentpart(p, rsendpartbuff, rsendpos, partpos)
      ELSE IF(sendholes(k) .LT. 0) THEN
        lsendpos=lsendpos+1
        CALL Insertsentpart(p, lsendpartbuff, lsendpos, partpos)
      END IF
    END DO
!
!
  END SUBROUTINE AddPartSendBuffers

  SUBROUTINE add_list_created_part(p, buffer,nb_ins)
    IMPLICIT NONE
    TYPE(particles), INTENT(INOUT):: p 
    TYPE(particle), ALLOCATABLE, INTENT(in) :: buffer(:)
    INTEGER, OPTIONAL:: nb_ins
    INTEGER:: i, nptotinit, parts_size_increase, nb_added

    nptotinit=p%Nploc+1
    if(present(nb_ins)) THEN
      nb_added=nb_ins
    ELSE
      nb_added=size(buffer,1)
    end if

    IF(nb_added .le. 0) RETURN ! No particles to add

    ! if there is not enough space in the parts simulation buffer, increase the parst size
    IF(p%Nploc + nb_added .gt. size(p%Z,1)) THEN
      parts_size_increase=Max(floor(0.1*size(p%Z,1)),nb_added)
      CALL change_parts_allocation(p, parts_size_increase)
    END IF
    
    DO i=1,nb_added
      CALL add_created_particle(p,buffer(i))
    END DO
    nb_added=p%Nploc-nptotinit+1
    if(p%is_field) then
      IF(allocated(p%addedlist)) then
        call change_array_size_int(p%addedlist,2)
      else
        allocate(p%addedlist(2))
      end if
      p%addedlist(size(p%addedlist)-1)=nptotinit
      p%addedlist(size(p%addedlist))=nb_added
    end if
  END SUBROUTINE add_list_created_part

  SUBROUTINE add_linked_created_part(p, linked_buffer, destroy, zerovelocity)

    IMPLICIT NONE
    TYPE(particles), INTENT(INOUT):: p 
    TYPE(linked_part_row), INTENT(in) :: linked_buffer
    LOGICAL:: destroy, zerovelocity
    TYPE(linked_part), POINTER:: part
    INTEGER:: i, nptotinit, parts_size_increase, nb_added

    nptotinit=p%Nploc+1
    nb_added=linked_buffer%n

    IF(nb_added .le. 0) RETURN ! No particles to add

    ! if there is not enough space in the parts simulation buffer, increase the parst size
    IF(p%Nploc + nb_added .gt. size(p%Z,1)) THEN
      parts_size_increase=Max(floor(0.1*size(p%Z,1)),nb_added)
      CALL change_parts_allocation(p, parts_size_increase)
    END IF
    
    part=>linked_buffer%start
    DO i=1,nb_added
      CALL add_created_particle(p,part%p)
      part=>part%next
    END DO
    nb_added=p%Nploc-nptotinit+1
    if(p%is_field) then
      IF(allocated(p%addedlist)) then
        call change_array_size_int(p%addedlist,2)
      else
        allocate(p%addedlist(2))
      end if
      p%addedlist(size(p%addedlist)-1)=nptotinit
      p%addedlist(size(p%addedlist))=nb_added
    end if
    if(zerovelocity)then
      p%UR(nptotinit:p%Nploc)=0
      p%UTHET(nptotinit:p%Nploc)=0
      p%UZ(nptotinit:p%Nploc)=0
    end if
    if (destroy) call destroy_linked_parts(linked_buffer%start)
    if (p%is_field) then
      ! we keep track of energy by removing the ionisation energy
      ! with conversion from electronvolt to joules
      loc_etot0=loc_etot0-sum(p%pot(nptotinit:p%Nploc)*elchar) 
    end if
  END SUBROUTINE add_linked_created_part

  !---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Copy the particles from the local parts variable to the left and right send buffers.
!
!> @param [in] lsendnbparts number of particles to send to the left neighbour (mpirank-1)
!> @param [in] rsendnbparts number of particles to send to the right neighbour (mpirank+1)
!> @param [in] sendholes array containing the indices of the particle leaving the local domain in ascending order. If the index is positive, the particle goes to the right neigbour, and to the left neighbour if the index is negative
!---------------------------------------------------------------------------

  SUBROUTINE add_created_particle(p,part)
    USE geometry
    TYPE(particles):: p
    TYPE(particle):: part
     p%Nploc=p%Nploc+1
      p%newindex=p%newindex+1
      ! add the data to the p structure
      CALL Insertincomingpart(p, part, p%Nploc)
      p%partindex(p%Nploc)=p%newindex
      ! calculate the new domain weight
      CALL dom_weight(p%Z(p%Nploc),p%R(p%Nploc),p%geomweight(p%Nploc,0))
      ! delete the particle if it is outside of the computational domain
      if( .not. is_inside(p,p%Nploc) ) then
        p%Nploc=p%Nploc-1
        p%newindex=p%newindex-1
        RETURN
      end if
      ! Calculate the geometric weight for the Poisson solver and the grid indices
      CALL geom_weight(p%Z(p%Nploc),p%R(p%Nploc),p%geomweight(p%Nploc,:))
      call p_calc_rzindex(p,p%Nploc)
    END SUBROUTINE add_created_particle

 function is_inside(p,id)
  Use basic, ONLY: rgrid,zgrid, nr, nz
  IMPLICIT NONE
   logical :: is_inside
   type(particles) :: p
   integer :: id   
   is_inside=.true.
   if(p%geomweight(id,0).le.0)then
      is_inside=.false.
      return
   end if
   if(p%R(id).ge.rgrid(nr) .or. p%R(id) .le. rgrid(0))then
    is_inside=.false.
    return
   end if
   if(p%Z(id).ge.zgrid(nz) .or. p%Z(id) .le. zgrid(0))then
    is_inside=.false.
    return
   end if
  end function is_inside

  SUBROUTINE calc_newparts_energy(p)
    USE basic, ONLY: phinorm, nlclassical
    type(particles)::p
    integer::i,n,nptotinit,nbadded, nptotend
    if(.not. p%is_field) return
    if( allocated(p%addedlist)) then
      n=size(p%addedlist)
      !write(*,*) n, "addedlist: ", p%addedlist
      Do i=1,n,2
        nptotinit=p%addedlist(i)
        nbadded=p%addedlist(i+1)
        p%nbadded=p%nbadded+nbadded
        nptotend=nptotinit+nbadded-1

        ! Potential energy
        loc_etot0=loc_etot0+p%q*p%weight*sum(p%pot(nptotinit:nptotend))*phinorm

        ! Kinetic energy
        IF(.not. nlclassical) THEN
            loc_etot0=loc_etot0+p%m*p%weight*vlight**2*sum(0.5*(p%Gamma(nptotinit:nptotend)+p%Gammaold(nptotinit:nptotend))-1)
        ELSE
            loc_etot0=loc_etot0+0.5*p%m*p%weight*vlight**2*sum(p%UR(nptotinit:nptotend)*p%URold(nptotinit:nptotend) &
            & +p%UZ(nptotinit:nptotend)*p%UZold(nptotinit:nptotend) &
            & +p%UTHET(nptotinit:nptotend)*p%UTHETold(nptotinit:nptotend))
        END IF
      end do
      deallocate(p%addedlist)
    end if

  end subroutine calc_newparts_energy

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Delete particle at given index removing its energy from the system
!
!> @param [in] index index of particle to be deleted
!---------------------------------------------------------------------------
  SUBROUTINE delete_part(p, index, replace)
  !! This will destroy particle at the given index
  USE constants, ONLY: vlight
  USE bsplines
  USE geometry
  USE basic, ONLY: phinorm, nlclassical
  TYPE(particles), INTENT(INOUT):: p
  INTEGER, INTENT(IN) :: index
  LOGICAL, OPTIONAL :: replace
  LOGICAL:: repl

  IF(present(replace)) THEN
    repl=replace
  ELSE
    repl=.true.
  END IF
  !Computes the potential at the particle position with phi_ext+phi_s
    IF(index .le. p%Nploc) THEN
      IF(p%is_field) THEN
        loc_etot0=loc_etot0-p%q*p%weight*(p%pot(index))*phinorm
        IF(.not. nlclassical) THEN
          loc_etot0=loc_etot0-p%m*p%weight*vlight**2*(p%Gamma(index)-1)
        ELSE
          loc_etot0=loc_etot0-0.5*p%m*p%weight*vlight**2*(p%UR(index)**2+p%UZ(index)**2+p%UTHET(index)**2)
        END IF
      END IF

      IF(repl) THEN
       ! We fill the gap
        CALL move_part(p, p%Nploc, index)
        p%partindex(p%Nploc)=-1
       ! Reduce the total number of simulated parts
        p%Nploc=p%Nploc-1
      END IF
    END IF
  END SUBROUTINE delete_part
!_______________________________________________________________________________
  SUBROUTINE loaduniformRZ(p, VR,VZ,VTHET)
    USE basic, ONLY: plasmadim, rnorm, temp, qsim, msim
    USE constants, ONLY: me, kb, elchar
    REAL(kind=db), INTENT(inout) ::VZ(:), VR(:), VTHET(:)
    TYPE(particles), INTENT(INOUT):: p

    CALL creat_parts(p, size(VR,1))

    p%Nploc=size(VR,1)
    p%Nptot=size(VR,1)
    p%q=sign(elchar,qsim)
    p%weight=msim/me
    p%m=me
    p%qmRatio=qsim/msim

  ! Initial distribution in z with normalisation
      CALL loduni(1,p%Z(1:p%Nploc))
      p%Z(1:p%Nploc)=(plasmadim(1)+(plasmadim(2)-plasmadim(1))*p%Z(1:p%Nploc))/rnorm
    ! Initial distribution in r with normalisation
      CALL lodlinr(2,p%R(1:p%Nploc),plasmadim(3),plasmadim(4))
      p%R(1:p%Nploc)=p%R(1:p%Nploc)/rnorm

    ! Initial velocities distribution
      CALL loadGaussianVelocities(p, VR, VZ, VTHET, temp)
  END SUBROUTINE loaduniformRZ
!_______________________________________________________________________________
  SUBROUTINE loadDavidson(p, VR,VZ,VTHET, lodr)
    USE constants, ONLY: me, kb, elchar
    USE basic, ONLY: nplasma, rnorm, plasmadim, distribtype, H0, P0, Rcurv, width, qsim, msim, &
    &                omegac, zgrid, nz, rnorm, n0, nblock, temp
    procedure(rloader)::lodr
    TYPE(particles), INTENT(INOUT):: p
    REAL(kind=db), INTENT(INOUT)::VZ(:), VR(:), VTHET(:)
    REAL(kind=db), DIMENSION(:), ALLOCATABLE::ra, rb, z
    REAL(kind=db) :: r0, deltar2, halfLz, Mirrorratio, Le, VOL
    INTEGER :: j, n, blockstart, blockend, addedpart, remainparts
    INTEGER, DIMENSION(:), ALLOCATABLE ::  blocksize

    CALL creat_parts(p, size(VR,1))

    p%Nploc=size(VR,1)
    p%Nptot=p%Nploc

      Allocate(ra(nblock),rb(nblock), z(0:nblock))
      !r0=(plasmadim(4)+plasmadim(3))/2
      r0=sqrt(4*H0/(me*omegac**2))
      halfLz=(zgrid(nz)+zgrid(0))/2
      MirrorRatio=(Rcurv-1)/(Rcurv+1)
      z(0)=plasmadim(1)
      DO n=1,nblock
      ! Compute limits in radius and load radii for each part
        Le=(plasmadim(2)-plasmadim(1))/nblock*(n-0.5)-halfLz*rnorm+plasmadim(1)
        z(n)=z(0)+n*(plasmadim(2)-plasmadim(1))/nblock
        deltar2=1-MirrorRatio*cos(2*pi*Le/width)
        rb(n)=r0/deltar2*sqrt(1-P0*abs(omegac)/2/H0*deltar2+sqrt(1-P0*abs(omegac)/H0*deltar2))
        ra(n)=r0/deltar2*sqrt(1-P0*abs(omegac)/2/H0*deltar2-sqrt(1-P0*abs(omegac)/H0*deltar2))
      END DO

      VOL=SUM(2*pi*MINVAL(ra)*(rb-ra)*(plasmadim(2)-plasmadim(1))/nblock)
      qsim=VOL*n0*elchar/nplasma
      msim=abs(qsim)/elchar*me
      p%weight=abs(qsim)/elchar
      p%m=me
      p%q=sign(elchar,qsim)
      p%qmRatio=p%q/p%m


      blockstart=1
      blockend=0
      ALLOCATE(blocksize(nblock))
      WRITE(*,*) "blocksize: ", size(blocksize), nblock
      DO n=1,nblock
        blocksize(n)=nplasma/VOL*2*pi*MINVAL(ra)*(rb(n)-ra(n))*(plasmadim(2)-plasmadim(1))/nblock
      END DO
      remainparts=p%Nploc-SUM(blocksize)
      addedpart=1
      n=nblock/2
      j=1
      DO WHILE(remainparts .GT. 0)
         blocksize(n)=blocksize(n)+addedpart
         remainparts=remainparts-addedpart
         n=n+j
         j=-1*(j+SIGN(1,j))
      END DO

      CALL loadPartSlices(p, lodr, ra, rb, z, blocksize)


    IF(distribtype .eq. 5) THEN
      CALL loadGaussianVelocities(p, VR, VZ, VTHET, temp)
      VZ=VZ/4
      VR=VR*8
      VTHET=VTHET*8
    ELSE
      Call loadDavidsonVelocities(p, VR, VZ, VTHET, H0, P0)
    END IF

  END SUBROUTINE loadDavidson


  SUBROUTINE loadDavidsonVelocities(p, VR,VZ,VTHET, H0, P0)
    USE constants, ONLY: me, kb, elchar
    USE basic, ONLY: rnorm, Rcurv, B0, width, vnorm, zgrid, nz
    TYPE(particles), INTENT(INOUT):: p
    REAL(kind=db), INTENT(INOUT)::VZ(:), VR(:), VTHET(:)
    REAL(kind=db), INTENT(IN):: H0, P0 
    REAL(kind=db) :: athetpos, rg, zg, halfLz, Mirrorratio, Pcomp, Acomp
    INTEGER :: i

    MirrorRatio=(Rcurv-1)/(Rcurv+1)
    halfLz=(zgrid(nz)+zgrid(0))/2

    ! Load velocities theta velocity
    ! Loading of r and z velocity is done in adapt_vinit to have
    ! access to parts%pot
      DO i=1,p%Nploc
    ! Interpolation for Magnetic potential
       rg=p%R(i)*rnorm
       zg=(p%Z(i)-halfLz)*rnorm

       Athetpos=0.5*B0*(rg - width/pi*MirrorRatio*bessi1(2*pi*rg/width)*COS(2*pi*zg/width))
       Pcomp=P0/rg/p%m
       Acomp=-p%qmRatio*Athetpos
       VTHET(i)=SIGN(MIN(abs(Pcomp+Acomp),sqrt(2*H0/p%m)),Pcomp+Acomp)
       !VTHET(i)=Pcomp+Acomp
      END DO
      VTHET=VTHET/vnorm
      VZ=0._db
      VR=0._db
      p%Davidson=.true.
      p%H0=H0
      p%P0=P0
  END SUBROUTINE loadDavidsonvelocities

  SUBROUTINE loadGaussianVelocities(p, VR,VZ,VTHET, temperature)
    USE basic, ONLY: vnorm
    USE constants, ONLY: kb
    REAL(kind=db), INTENT(inout) ::VZ(:), VR(:), VTHET(:)
    TYPE(particles), INTENT(INOUT):: p
    REAL(kind=db), INTENT(IN):: temperature
    REAL(kind=db):: vth

    ! Initial velocities distribution
      vth=sqrt(2.0/3.0*kb*temperature/p%m)/vnorm        !thermal velocity
      CALL lodgaus(3,VZ)
      CALL lodgaus(5,VR)
      CALL lodgaus(7,VTHET)
      VZ=VZ*vth
      VR=VR*vth
      VTHET=VTHET*vth
      p%temperature=temperature
      p%Davidson=.false.
  END SUBROUTINE loadGaussianVelocities

  SUBROUTINE loadFlatTopVelocities(p, VR,VZ,VTHET, meanv, spanv)
    USE basic, ONLY: vnorm
    USE constants, ONLY: kb
    REAL(kind=db), INTENT(inout) ::VZ(:), VR(:), VTHET(:)
    TYPE(particles), INTENT(INOUT):: p
    REAL(kind=db), INTENT(INOUT):: meanv(3), spanv(3) 

    ! Initial velocities distribution
      meanv=meanv/vnorm        !thermal velocity
      spanv=spanv/vnorm
      CALL loduni(3,VZ)
      CALL loduni(5,VR)
      CALL loduni(7,VTHET)
      VR=(VR*2-1)*spanv(1)+meanv(1)
      VTHET=(VTHET*2-1)*spanv(2)+meanv(2)
      VZ=(VZ*2-1)*spanv(3)+meanv(3)
      p%Davidson=.false.
  END SUBROUTINE loadFlatTopVelocities


  SUBROUTINE loadPartslices(p, lodr, ra, rb, z, blocksize)
    USE basic, ONLY: rnorm
    TYPE(particles), INTENT(INOUT):: p
    REAL(kind=db), INTENT(IN)::ra(:), rb(:), z(0:)
    INTEGER, DIMENSION(:), INTENT(IN) ::  blocksize
    procedure(rloader)::lodr
    INTEGER :: n, blockstart, blockend, nblock
      nblock=size(blocksize,1)
      blockstart=1
      blockend=0
      DO n=1,nblock
        blockstart=blockend+1
        blockend=MIN(blockstart+blocksize(n)-1,p%Nploc)
    ! Initial distribution in z with normalisation between magnetic mirrors
        CALL loduni(1, p%Z(blockstart:blockend))
        p%Z(blockstart:blockend)= (z(n-1)+p%Z(blockstart:blockend)*(z(n)-z(n-1)))/rnorm
        CALL lodr(2, p%R(blockstart:blockend), ra(n), rb(n))
        p%R(blockstart:blockend)=p%R(blockstart:blockend)/rnorm
      END DO
  END SUBROUTINE loadPartslices


  !---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Read a particle file format to load a simulated specie in the simulation
!
!---------------------------------------------------------------------------
  SUBROUTINE read_part_file(p, partfilename, VR, VZ, VTHET)
    USE basic, ONLY: lu_partfile, rnorm, vnorm
    implicit None
    TYPE(particles), INTENT(INOUT):: p
    REAL(kind=db), DIMENSION(:), ALLOCATABLE, INTENT(INOUT)::VR, VZ, VTHET
    CHARACTER(len=*)::partfilename
    INTEGER:: nblock = 0
    REAL(kind=db), Dimension(:), ALLOCATABLE:: ra, rb, z
    INTEGER, Dimension(:), ALLOCATABLE:: npartsslice
    INTEGER:: velocitytype=1  !< 1) gaussian with temp  2) Davidson with H0, P0  
    INTEGER:: radialtype=1    !< 1) 1/R  2) uniform  3) 1/R^2 4) gauss
    INTEGER:: npartsalloc     !< initial size of particles arrays
    INTEGER:: iiee_id         !< index of species to add particles to for IIEE
    INTEGER:: neuttype_id     !< index of neutral gas producing ions
    REAL(kind=db):: mass=me
    REAL(kind=db):: charge=-elchar
    REAL(kind=db):: weight=1.0
    REAL(kind=db):: qmratioscale
    REAL(kind=db):: meanv(3)          !< mean velocity in each direction for velocitytype 3
    REAL(kind=db):: spanv(3)          !< pos/neg extent of velocity in each direction for velocitytype 3
    CHARACTER(len=256) :: header=' '  !< header of csv file
    REAL(kind=db):: H0=3.2e-14        !< Total energy
    REAL(kind=db):: P0=8.66e-25       !< Canonical angular momentum
    REAL(kind=db):: temperature=10000 !< temperature in kelvins
    real(kind=db):: n0                !< density factor
    LOGICAL :: is_test                !< Defines if particle are saved on ittracer or not 
    LOGICAL :: is_field               !< Defines if particle contributes to Poisson solver
    LOGICAL :: calc_moments           !< Defines if moments matrix must be calculated each it2d
    CHARACTER(len=16) :: partformat = 'slice'
    INTEGER:: i, ierr, openerr

    NAMELIST /partsload/ nblock, mass, charge, weight, npartsalloc, velocitytype, & 
             & radialtype, temperature, H0, P0, is_test, n0, partformat, meanv, spanv, &
             & calc_moments, qmratioscale, is_field, iiee_id, neuttype_id
    
    ! Set defaults
    qmratioscale=1.0
    weight=1.0
    meanv=0
    spanv=0
    mass=me
    charge=-elchar
    calc_moments=.false.
    is_test=.false.
    is_field=.true.
    iiee_id = -1
    neuttype_id=1

    ! Open the paticle file
    OPEN(UNIT=lu_partfile,FILE=trim(partfilename),ACTION='READ',IOSTAT=openerr)
    header=' '
    IF(openerr .ne. 0) THEN
      CLOSE(unit=lu_partfile)
      RETURN
    END IF
    READ(lu_partfile,partsload)
    
    IF(mpirank .eq.0) THEN
      WRITE(*,'(a,a)')"reading partfile: ", trim(partfilename) 
      WRITE(*,partsload)
    END IF
   
    ! The plasma cloud is defined as a set of slices 
    IF(trim(partformat).eq.'slice') THEN
      IF( nblock .ge. 1) THEN
        ALLOCATE(z(0:nblock),ra(nblock),rb(nblock), npartsslice(nblock))
        DO WHILE(header(1:8) .ne. '//slices')
          READ(lu_partfile,'(a)') header
        END DO
        DO i=1,nblock
          READ(lu_partfile,*) z(i-1),ra(i),rb(i),npartsslice(i)
        END DO
        READ(lu_partfile,*) z(nblock)

        CALL creat_parts(p,max(npartsalloc,sum(npartsslice)))
        p%Nploc=sum(npartsslice)
        p%Nptot=p%Nploc
        IF( allocated(VR) ) THEN
          DEALLOCATE(VR,VZ,VTHET)
        end if
        if(.not. allocated(VR)) THEN
          ALLOCATE(VR(p%Nploc))
          ALLOCATE(VZ(p%Nploc))
          ALLOCATE(VTHET(p%Nploc))
        END IF
        
        p%m=mass
        p%q=charge
        p%weight=weight
        p%qmRatio=charge/mass*qmratioscale
        p%is_test=is_test
        p%is_field=is_field
        p%calc_moments=calc_moments
        p%Newindex=sum(npartsslice)
        p%iiee_id = iiee_id
        p%neuttype_id = neuttype_id

        SELECT CASE(radialtype)
          CASE(1) ! 1/R distribution in R
            CALL loadPartslices(p, lodunir, ra, rb, z, npartsslice)
          CASE(2) ! flat top distribution in R
            CALL loadPartslices(p, lodlinr, ra, rb, z, npartsslice)
          CASE(3) ! 1/R^2 distribution in R
            CALL loadPartslices(p, lodinvr, ra, rb, z, npartsslice)
          CASE(4) ! gaussian distribution in R
            CALL loadPartslices(p, lodgausr, ra, rb, z, npartsslice)
          CASE DEFAULT
            IF (mpirank .eq. 0) WRITE(*,*) "Unknown type of radial distribution:", radialtype
            CALL MPI_Abort(MPI_COMM_WORLD, -1, ierr)
        END SELECT

        SELECT CASE(velocitytype)
        CASE(1) ! Gaussian with temperature
          CALL loadGaussianVelocities(p, VR, VZ, VTHET, temperature)
        CASE(2) ! Davidson magnetic mirror high wr equilibrium
          CALL loadDavidsonVelocities(p, VR, VZ, VTHET, H0, P0)
        CASE(3) ! flat top velocity
          CALL loadFlatTopVelocities(p, VR, VZ, VTHET, meanv, spanv)
        CASE DEFAULT
          IF (mpirank .eq. 0) WRITE(*,*) "Unknown type of velocity distribution:", velocitytype
          CALL MPI_Abort(MPI_COMM_WORLD, -1, ierr)
        END SELECT

      END IF

    END IF


    ! The plasma cloud is defined as a set individual particles
    IF( trim(partformat) .eq. 'parts' ) THEN
      IF( nblock .ge. 1) THEN
        !Allocate necessary memory
        CALL creat_parts(p,max(npartsalloc,nblock))
        IF( allocated(VR) ) THEN
          DEALLOCATE(VR,VZ,VTHET)
        end if
        if(.not. allocated(VR)) THEN
          ALLOCATE(VR(nblock))
          ALLOCATE(VZ(nblock))
          ALLOCATE(VTHET(nblock))
        END IF

        ! Read the particles from the file
        DO WHILE(header(1:8) .ne. '//parts')
          READ(lu_partfile,'(a)') header
        END DO
        DO i=1,nblock
          READ(lu_partfile,*) p%R(i),p%THET(i),p%Z(i), VR(i), VTHET(i), VZ(i)
        END DO

        p%Nploc=nblock
        p%Nptot=p%Nploc
        p%m=mass
        p%q=charge
        p%Newindex=nblock
        p%weight=weight
        p%qmRatio=charge/mass*qmratioscale
        p%is_test=is_test
        p%is_field=is_field
        p%calc_moments=calc_moments
        p%iiee_id = iiee_id
        p%neuttype_id = neuttype_id
        
        !normalizations
        p%r=p%r/rnorm
        p%z=p%z/rnorm
        VR=VR/vnorm
        VTHET=VTHET/vnorm
        VZ=VZ/vnorm
      END IF

    END IF
    CLOSE(unit=lu_partfile)
    


  END SUBROUTINE

!---------------------------------------------------------------------------
!> @author
!> Guillaume Le Bars EPFL/SPC
!
! DESCRIPTION:
!>
!> @brief Increase the number of macroparticles by separating each previous macroparticles into 
!> samplefactor new macroparticles of equally divided weight. The new sub particles are distributed 
!> uniformly in space to maintain the density and other moments.
!
!> @param [in] samplefactor multiplicator of the number of macroparticles.
!> @param [in] p particles type to increase.
!---------------------------------------------------------------------------

  SUBROUTINE upsample(p, samplefactor)
    USE basic, ONLY : nplasma, dr, dz
    INTEGER, INTENT(IN) ::samplefactor
    TYPE(particles), INTENT(INOUT):: p
    INTEGER:: i, j, currentindex
    REAL(kind=db), DIMENSION(p%Nploc) :: spreaddir ! random direction for the spread of each initial macro particle
    REAL(kind=db) :: dir ! Direction in which the particle is moved
    REAL(kind=db) :: dl  ! Particle displacement used for
    ! Load and scale the direction angle for spreading the new particles
    CALL loduni(2, spreaddir)
    spreaddir=spreaddir*2*pi/samplefactor
    dl=min(dz,minval(dr,1,dr.GT.0))/100
    DO i=1,p%Nploc
      DO j=1,samplefactor-1
        currentindex=p%Nploc+(i-1)*(samplefactor-1)+j
        CALL move_part(p,i,currentindex)
        p%partindex(currentindex)=currentindex
        dir = spreaddir(i)+2*pi*j/samplefactor
        p%R(currentindex)=p%R(currentindex) + dl*cos(dir)
        p%Z(currentindex)=p%Z(currentindex) + dl*sin(dir)
      END DO
      p%partindex(i)=i
      p%R(i)=p%R(i) + dl*cos(spreaddir(i))
      p%Z(i)=p%Z(i) + dl*sin(spreaddir(i))
    END DO
    nplasma=nplasma*samplefactor
    p%weight=p%weight/samplefactor
    p%Nploc=p%Nploc*samplefactor
    p%Nptot=p%Nptot*samplefactor
  END SUBROUTINE upsample


! Taken from https://rosettacode.org/wiki/Sorting_algorithms/Radix_sort#Fortran
!	No Copyright is exerted due to considerable prior art in the Public Domain.
!       This Fortran version by Peter Kelly ~ peter.kelly@acm.org
!
!	Permission is hereby granted, free of charge, to any person obtaining
!	a copy of this software and associated documentation files (the
!	"Software"), to deal in the Software without restriction, including
!	without limitation the rights to use, copy, modify, merge, publish,
!	distribute, sublicense, and/or sell copies of the Software, and to
!	permit persons to whom the Software is furnished to do so, subject to
!	the following conditions:
!	The above copyright notice and this permission notice shall be
!	included in all copies or substantial portions of the Software.
!	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!	CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!	TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!          
!     Implementation of a classic Radix Sort LSD style :)
  SUBROUTINE LSDRADIXSORT(A , N)
    IMPLICIT NONE
!
! Dummy arguments
!
    INTEGER  ::  N
    INTEGER , target, DIMENSION(0:N - 1)  ::  A           ! All arrays based off zero, one day I'll fix it
    INTENT (IN) N
    INTENT (INOUT) A
!
! Local variables
! 
    INTEGER , DIMENSION(0:9)  ::  counts
    INTEGER  ::  digitplace
    INTEGER  ::  i
    INTEGER  ::  j
    INTEGER  ::  largestnum
    INTEGER, DIMENSION(0:N - 1)  ::  results 
! 
    digitplace = 1                                        ! Count of the keys
    largestnum = MAXVAL(A)

    DO WHILE ( (largestnum/digitplace)>0 )
       counts = 0                                         ! Init the count array
      DO i = 0 , N - 1 , 1
          J = (A(i)/digitplace)
          J = MODULO(j , 10) 
          counts(j) = counts(j) + 1
      END DO

!  Change count(i) so that count(i) now contains actual position of this digit in result()
!  Working similar to the counting sort algorithm
       DO i = 1 , 9 , 1
          counts(i) = counts(i) + counts(i - 1)       ! Build up the prefix sum
       END DO
!
       DO i = N - 1 , 0 , -1                          ! Move from left to right
          j = (A(i)/digitplace)
          j = MODULO(j, 10)
          results(counts(j) - 1) = A(i)               ! Need to subtract one as we are zero based but prefix sum is 1 based
          counts(j) = counts(j) - 1
       END DO
!
       DO i = 0 , N - 1 , 1                           ! Copy the semi-sorted data into the input
         A(i) = results(i)
       END DO
!
       digitplace = digitplace*10
    END DO                                             ! While loop
    RETURN
    END SUBROUTINE LSDRADIXSORT

END MODULE beam
