SUBROUTINE chkrst(flag)
!
!   Process checkpoint/restart file
!
  USE basic
  USE futils
  USE beam
  USE fields
  USE constants, ONLY: elchar, me
  Use psupply
  IMPLICIT NONE
  INTEGER, INTENT(in) :: flag
  INTEGER :: remainingparts
  REAL(kind=db):: old_msim, old_qsim, old_n0
  INTEGER:: partsrank, partsdim(1), i, err
  REAL(kind=db), ALLOCATABLE:: charges(:), weights(:), masses(:) 
  CHARACTER(len=64):: group
  CHARACTER(len=2):: specieindex
  real(kind=db):: old_rnorm, old_tnorm, qmratioscale
  real(kind=db):: prev_bias
  INTEGER:: logical_val, id
  ! Only process 0 should save on file
  
 
  !
  !  Local vars and arrays
  !________________________________________________________________________________
  !
    SELECT CASE(flag)
  !________________________________________________________________________________
  !                   1.  Open and read restart file
  !
    CASE(0)
      CALL openf(rstfile, fidrst,'r',real_prec='d')
      CALL getatt(fidrst, '/Basic', 'cstep', cstep)
      CALL getatt(fidrst, '/Basic', 'time', time)
      CALL getatt(fidrst, '/Basic', 'n0',   old_n0)

      IF(isgroup(fidrst,'/Basic/norm')) THEN
        CALL getatt(fidrst, '/Basic/norm', 'rnorm',   old_rnorm)
        CALL getatt(fidrst, '/Basic/norm', 'tnorm',   old_tnorm)
      else
        old_rnorm=rnorm
        old_tnorm=tnorm
      end if
      
      IF(isdataset(fidrst,'/Parts/charges')) THEN
      ! If we have multiple saved species we need to load all of them
        CALL getatt(fidrst,'/Parts','nbspecies',nbspecies)
        nbspecies=nbspecies
        ALLOCATE(charges(nbspecies),masses(nbspecies),weights(nbspecies))
        ALLOCATE(partslist(nbspecies+nbaddtestspecies))
        CALL getarr(fidrst, '/Parts/charges', charges)
        CALL getarr(fidrst, '/Parts/masses',  masses)
        CALL getarr(fidrst, '/Parts/weights', weights)
        weights(1)=weights(1)/old_n0*n0
      ELSE
      ! If we have an old restart file, we load only the electrons
        CALL getatt(fidrst, '/Basic', 'msim', old_msim)
        CALL getatt(fidrst, '/Basic', 'qsim', old_qsim)
        qsim=old_qsim/old_n0*n0
        msim=old_msim/old_n0*n0
        nbspecies=1
        ALLOCATE(charges(nbspecies),masses(nbspecies),weights(nbspecies))
        ALLOCATE(partslist(nbspecies))
        charges(1)=sign(elchar,qsim)
        weights(1)=msim/me
        masses(1)=me
      END IF
      
      if(newres) then
        weights=weights*weights_scale
      end if

      CALL getatt(fidrst, '/var0d', 'etot0', loc_etot0)
      etot0=loc_etot0
      if(n0.ne.old_n0) cstep=0
      loc_etot0=0
      CALL getatt(fidrst, '/var0d', 'epot', epot)
      CALL getatt(fidrst, '/var0d', 'ekin', ekin)
      CALL getatt(fidrst, '/var0d', 'etot', etot)

      CALL getatt(fidrst, '/Parts', 'nplasma', nplasma)
      CALL getatt(fidrst, '/Parts', 'remainingparts', remainingparts)
      if(remainingparts .gt. 0) Then
        CALL getdims(fidrst, '/Parts/Z', partsrank, partsdim)
      else
        partsdim=0
      end if
      IF (samplefactor .gt. 1 ) THEN ! We increase the number of macro particles
        CALL creat_parts(partslist(1),max(remainingparts*samplefactor,partsdim(1)))
      ELSE
        CALL creat_parts(partslist(1),partsdim(1))
      END IF
      partslist(1)%q=charges(1)
      partslist(1)%m=masses(1)
      partslist(1)%weight=weights(1)
      partslist(1)%qmRatio=charges(1)/masses(1)
      err=0
      CALL getatt(fidrst, 'Parts', 'qmratioscale', qmratioscale, err)
      if(err .ge.0) partslist(1)%qmRatio=partslist(1)%qmRatio*qmratioscale

      if(remainingparts .gt. 0) then
        CALL getarr(fidrst, '/Parts/Z', partslist(1)%Z)
        CALL getarr(fidrst, '/Parts/R', partslist(1)%R)
        ! Renormalize R and Z
        IF(isgroup(fidrst,'/Basic/norm')) THEN
          partslist(1)%R=partslist(1)%R*old_rnorm/rnorm
          partslist(1)%Z=partslist(1)%Z*old_rnorm/rnorm
        ELSE
          partslist(1)%R=partslist(1)%R*sqrt(n0/old_n0)
          partslist(1)%Z=partslist(1)%Z*sqrt(n0/old_n0)
        END IF
        CALL getarr(fidrst, '/Parts/THET', partslist(1)%THET)
        CALL getarr(fidrst, '/Parts/UZ', partslist(1)%UZ)
        CALL getarr(fidrst, '/Parts/UR', partslist(1)%UR)
        CALL getarr(fidrst, '/Parts/UTHET', partslist(1)%UTHET)
        CALL getarr(fidrst, '/Parts/UZ', partslist(1)%UZold)
        CALL getarr(fidrst, '/Parts/UR', partslist(1)%URold)
        CALL getarr(fidrst, '/Parts/UTHET', partslist(1)%UTHETold)
        CALL getarr(fidrst, '/Parts/Zindex', partslist(1)%Zindex)
        CALL getarr(fidrst, '/Parts/Rindex', partslist(1)%Rindex)
        CALL getarr(fidrst, '/Parts/partindex', partslist(1)%partindex)
        IF(isdataset(fidrst,'/Parts/fluidur')) THEN
          CALL getarr(fidrst, '/Parts/GAMMA', partslist(1)%Gamma)
        END IF
      end if
      partslist(1)%Nploc=remainingparts
      partslist(1)%Nptot=partslist(1)%Nploc
      partslist(1)%Newindex=maxval(partslist(1)%partindex)
      WRITE(*,*) "Read ", partslist(1)%Nploc, " particles out of ", remainingparts
      
      IF(nbspecies .gt. 1) THEN
        DO i=2,nbspecies
          WRITE(group,'(a,i2)')'/Parts/',i
          WRITE(specieindex,'(i2)') i
          partsdim=0
          CALL getatt(fidrst, trim(group), 'remainingparts', remainingparts)
          if(remainingparts .gt. 0) Then
            CALL getdims(fidrst, trim(group) // '/Z', partsrank, partsdim)
          else
            partsdim=0
          end if
          IF(partsdim(1).gt.remainingparts) THEN
            CALL creat_parts(partslist(i),partsdim(1))
            partslist(i)%Nploc=remainingparts
          ELSE
            CALL creat_parts(partslist(i),max(10,remainingparts))
          ENDIF
          partslist(i)%q=charges(i)
          partslist(i)%m=masses(i)
          partslist(i)%weight=weights(i)
          partslist(i)%qmRatio=charges(i)/masses(i)
          err=0
          CALL getatt(fidrst, trim(group), 'qmratioscale', qmratioscale, err)
          if(err .ge.0) partslist(i)%qmRatio=partslist(i)%qmRatio*qmratioscale
          partslist(i)%Nptot=remainingparts
          partslist(i)%Nploc=remainingparts
          partslist(i)%is_test =.false.
          partslist(i)%is_field =.false.
          partslist(i)%calc_moments =.false.
          err=0    
          CALL getatt(fidrst, trim(group), 'is_test', logical_val,err)
          if(err .ge.0)then
            if(logical_val.gt.0) partslist(i)%is_test =.true.
          end if
          err=0
          CALL getatt(fidrst, trim(group), 'is_field', logical_val,err)
          if(err .ge.0)then
            if(logical_val.gt.0) partslist(i)%is_field =.true.
          end if
          err=0
          CALL getatt(fidrst, trim(group), 'calc_moments', logical_val,err)
          if(err .ge.0)then
            if(logical_val.gt.0) partslist(i)%calc_moments =.true.
          end if
          ! ---------------------------------------------------------------------------------
          ! IIEE PARAMETERS 
          CALL getatt(fidrst, trim(group), 'zero_vel', logical_val,err)
          if(err .ge.0)then
            if(logical_val.gt.0) partslist(i)%zero_vel =.true.
          end if
          CALL getatt(fidrst, trim(group), 'iiee_id',id, err)
          if (err .ge. 0) then 
             partslist(i)%iiee_id = id
          end if
          CALL getatt(fidrst, trim(group), 'neuttype_id',id, err)
          if (err .ge. 0) then
             partslist(i)%neuttype_id = id
          end if

          CALL getatt(fidrst, trim(group), 'material_id',id, err)
          if (err .ge. 0) then
             partslist(i)%material_id = id
          end if

          ! END IIEE PARAMETERS 
          ! ---------------------------------------------------------------------------------
          IF(partslist(i)%Nptot .gt. 0) THEN
            CALL getarr(fidrst, trim(group) // '/Z',         partslist(i)%Z)
            CALL getarr(fidrst, trim(group) // '/R',         partslist(i)%R)
            CALL getarr(fidrst, trim(group) // '/THET',      partslist(i)%THET)
            CALL getarr(fidrst, trim(group) // '/UZ',        partslist(i)%UZ)
            CALL getarr(fidrst, trim(group) // '/UR',        partslist(i)%UR)
            CALL getarr(fidrst, trim(group) // '/UTHET',     partslist(i)%UTHET)
            CALL getarr(fidrst, trim(group) // '/UZ',        partslist(i)%UZold)
            CALL getarr(fidrst, trim(group) // '/UR',        partslist(i)%URold)
            CALL getarr(fidrst, trim(group) // '/UTHET',     partslist(i)%UTHETold)
            CALL getarr(fidrst, trim(group) // '/GAMMA',     partslist(i)%Gamma)
            CALL getarr(fidrst, trim(group) // '/Zindex',    partslist(i)%Zindex)
            CALL getarr(fidrst, trim(group) // '/Rindex',    partslist(i)%Rindex)
            CALL getarr(fidrst, trim(group) // '/partindex', partslist(i)%partindex)
            IF(isgroup(fidrst,'/Basic/norm')) THEN
              partslist(i)%R=partslist(i)%R*old_rnorm/rnorm
              partslist(i)%Z=partslist(i)%Z*old_rnorm/rnorm
            ELSE
              partslist(i)%R=partslist(i)%R*sqrt(n0/old_n0)
              partslist(i)%Z=partslist(i)%Z*sqrt(n0/old_n0)
            END IF
            partslist(i)%Newindex=maxval(partslist(i)%partindex)
          END IF
        END DO
      END IF

      IF(isgroup(fidrst,'/psupply')) THEN
          call getatt(fidrst,'/psupply', 'active', logical_val)
          if(logical_val .gt. 0) then
            call getatt(fidrst,'/psupply', 'bias', prev_bias)
            the_ps%active=.true.
            the_ps%bias=prev_bias/phinorm
          else
            the_ps%active=.false.
          end if
      end if


      CALL closef(fidrst)
      IF(samplefactor .gt. 1) THEN ! We increase the number of macro particles
        CALL upsample(partslist(1), samplefactor)
      END IF
      WRITE(*,'(3x,a)') "Reading from restart file "//TRIM(rstfile)//" completed!"
  !________________________________________________________________________________
  !                   2.  Create and write to restart file (DP reals)
  !
    CASE(1)
       
      IF( .NOT. nlsave ) RETURN
      CALL mv2bk(rstfile)
      CALL creatf(rstfile, fidrst, real_prec='d', desc='Restart file')
      CALL creatg(fidrst, '/Basic', 'Basic data')
      CALL attach(fidrst, '/Basic', 'cstep', cstep)
      CALL attach(fidrst, '/Basic', 'time', time)
      CALL attach(fidrst, '/Basic', 'jobnum', jobnum)
      CALL attach(fidrst, '/Basic', 'qsim', partslist(1)%q*partslist(1)%weight)
      CALL attach(fidrst, '/Basic', 'msim', partslist(1)%m*partslist(1)%weight)
      CALL attach(fidrst, '/Basic', 'n0',   n0)
      CALL creatg(fidrst, '/Basic/norm', 'Normalisation quantities')
      CALL attach(fidrst, '/Basic/norm', 'rnorm',   1.0)
      CALL attach(fidrst, '/Basic/norm', 'bnorm',   bnorm)
      CALL attach(fidrst, '/Basic/norm', 'enorm',   enorm)
      CALL attach(fidrst, '/Basic/norm', 'tnorm',   tnorm)
      CALL attach(fidrst, '/Basic/norm', 'phinorm', phinorm)
    !
    !  0D variables
    !
      CALL creatg(fidrst, '/var0d', '0D variables')
      CALL attach(fidrst, '/var0d','etot0', etot0)
      CALL attach(fidrst, '/var0d','epot', epot)
      CALL attach(fidrst, '/var0d','ekin', ekin)
      CALL attach(fidrst, '/var0d','etot', etot)

    !
    !  Parts
    !
      CALL creatg(fidrst, '/Parts', 'Particles data')
      CALL attach(fidrst, '/Parts', 'nplasma', nplasma)
      nbspecies=size(partslist,1)
      CALL attach(fidrst, '/Parts', 'nbspecies', nbspecies)
      ALLOCATE(charges(nbspecies),masses(nbspecies),weights(nbspecies))

      DO i=1,nbspecies
        charges(i) = partslist(i)%q
        masses(i)  = partslist(i)%m
        weights(i) = partslist(i)%weight
      END DO
      CALL putarr(fidrst, '/Parts/charges', charges)
      CALL putarr(fidrst, '/Parts/masses',  masses )
      CALL putarr(fidrst, '/Parts/weights', weights)
      IF(mpisize .gt. 1) THEN
        remainingparts=sum(Nplocs_all)
      ELSE 
        remainingparts=partslist(1)%Nploc
      END IF
      CALL attach(fidrst, '/Parts', 'remainingparts', remainingparts)
      CALL attach(fidrst, '/Parts', 'qmratioscale',partslist(1)%qmRatio/(partslist(1)%q/partslist(1)%m))
      CALL putarr(fidrst, '/Parts/Z', partslist(1)%Z*rnorm)
      CALL putarr(fidrst, '/Parts/R', partslist(1)%R*rnorm)
      CALL putarr(fidrst, '/Parts/THET', partslist(1)%THET)
      CALL putarr(fidrst, '/Parts/UZ', partslist(1)%UZ)
      CALL putarr(fidrst, '/Parts/UR', partslist(1)%UR)
      CALL putarr(fidrst, '/Parts/UTHET', partslist(1)%UTHET)
      CALL putarr(fidrst, '/Parts/GAMMA', partslist(1)%Gamma)
      CALL putarr(fidrst, '/Parts/Zindex', partslist(1)%Zindex)
      CALL putarr(fidrst, '/Parts/Rindex', partslist(1)%Rindex)
      CALL putarr(fidrst, '/Parts/partindex', partslist(1)%partindex)
      CALL putarr(fidrst, '/Parts/fluidur', partslist(1)%moments(2,:))
      CALL putarr(fidrst, '/Parts/fluiduthet', partslist(1)%moments(3,:))
      CALL putarr(fidrst, '/Parts/fluiduz', partslist(1)%moments(4,:))
      partslist(1)%is_field=.true.
      partslist(1)%is_test=.false.
      partslist(1)%calc_moments=.true.
      CALL attach(fidrst, '/Parts', 'is_field', 1)
      CALL attach(fidrst, '/Parts', 'calc_moments', 1)
      CALL attach(fidrst, '/Parts', 'is_test', 0)
      

      IF(nbspecies .gt. 1) THEN
        DO i=2,nbspecies
          WRITE(group,'(a,i2)')'/Parts/',i
          WRITE(specieindex,'(i2)') i
          CALL creatg(fidrst, trim(group), 'Particles ' // specieindex// ' data')
          CALL attach(fidrst, trim(group), 'qmratioscale', partslist(i)%qmRatio/(partslist(i)%q/partslist(i)%m))
          CALL attach(fidrst, trim(group), 'remainingparts', partslist(i)%Nptot)
          if(partslist(i)%is_test)then
            CALL attach(fidrst, trim(group), 'is_test', 1)
          else
            CALL attach(fidrst, trim(group), 'is_test', 0)
          end if
          if(partslist(i)%is_field)then
            CALL attach(fidrst, trim(group), 'is_field', 1)
          else
            CALL attach(fidrst, trim(group), 'is_field', 0)
          end if
          if(partslist(i)%calc_moments)then
            CALL attach(fidrst, trim(group), 'calc_moments', 1)
          else
            CALL attach(fidrst, trim(group), 'calc_moments', 0)
          end if
          !-------------------------------------------------------------------------
          ! IIEE PARAMETERS 
          if(partslist(i)%zero_vel)then
            CALL attach(fidrst, trim(group), 'zero_vel', 1)
          else
            CALL attach(fidrst, trim(group), 'zero_vel', 0)
          end if
          CALL attach(fidrst, trim(group), 'material_id', partslist(i)%iiee_id)
          CALL attach(fidrst, trim(group), 'neuttype_id', partslist(i)%neuttype_id)
          CALL attach(fidrst, trim(group), 'material_id', partslist(i)%material_id)
          ! END OF IIEE PARAMETERS
          ! ------------------------------------------------------------------------
          IF(partslist(i)%Nptot .gt. 0) THEN
            CALL putarr(fidrst, trim(group) // '/Z',         partslist(i)%Z(1:partslist(i)%Nptot)*rnorm)
            CALL putarr(fidrst, trim(group) // '/R',         partslist(i)%R(1:partslist(i)%Nptot)*rnorm)
            CALL putarr(fidrst, trim(group) // '/THET',      partslist(i)%THET(1:partslist(i)%Nptot))
            CALL putarr(fidrst, trim(group) // '/UZ',        partslist(i)%UZ(1:partslist(i)%Nptot))
            CALL putarr(fidrst, trim(group) // '/UR',        partslist(i)%UR(1:partslist(i)%Nptot))
            CALL putarr(fidrst, trim(group) // '/UTHET',     partslist(i)%UTHET(1:partslist(i)%Nptot))
            CALL putarr(fidrst, trim(group) // '/GAMMA',     partslist(i)%Gamma(1:partslist(i)%Nptot))
            CALL putarr(fidrst, trim(group) // '/Zindex',    partslist(i)%Zindex(1:partslist(i)%Nptot))
            CALL putarr(fidrst, trim(group) // '/Rindex',    partslist(i)%Rindex(1:partslist(i)%Nptot))
            CALL putarr(fidrst, trim(group) // '/partindex', partslist(i)%partindex(1:partslist(i)%Nptot))
          END IF
        END DO
      END IF
    !
    !  Fields
    !

    !
    ! Power supply status
    !
      CALL creatg(fidrst, '/psupply', 'Power supply status and values')
      if(the_ps%active) then
        CALL attach(fidrst, '/psupply','active', 1)
      else
        CALL attach(fidrst, '/psupply','active', 0)
      end if
      CALL attach(fidrst, '/psupply','bias', the_ps%bias*phinorm)


      CALL closef(fidrst)
      WRITE(*,'(3x,a)') "Writing to restart file "//TRIM(rstfile)//" completed!"
  !
    END SELECT
!
 END SUBROUTINE chkrst
