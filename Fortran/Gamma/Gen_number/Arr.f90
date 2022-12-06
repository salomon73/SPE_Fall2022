PROGRAM init_table

        implicit none
        integer(kind = 4) :: nppts 
        real(kind = 8), dimension(10) :: tableau
        integer(kind = 4) :: ii
        real(kind = 8) :: Emin, Emax

        Emin = 0.0
        Emax = 25.0
        nppts = 10

        do ii = 1,size(tableau)
             tableau(ii) = (Emin + (ii-1)*(Emax-Emin)/(nppts-1) ) 
             print*, tableau(ii)
        end do 
END PROGRAM init_table

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

