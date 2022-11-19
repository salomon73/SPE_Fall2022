MODULE poisson_unif

USE factorial        
CONTAINS


INTEGER FUNCTION gen_elec(lambda, kmax)


    REAL :: lambda !< Lambda parameter for Poisson distribution 
    REAL :: nb_alea!< random number unif. generated in [0,1]
    INTEGER :: kmax           !< max number possible from Poisson
    REAL :: CumulPoisson  !< Flag to ensure CDF ~ 1
    INTEGER :: i, ii                 !< loop indices   
    REAL, DIMENSION(kmax) :: vect, SumPart !< terms, partial sums for CDF

    !> Compute probabilities for each int. value and CDF values
    DO i = 1,kmax
       vect(i)    = exp(-lambda)*lambda**(i-1)/factorial_fun(i-1);
       SumPart(i) = sum(vect(1:i));
    END DO

    !> CDF expected to sum to 1 
    CumulPoisson = sum(vect)

    !> Generate poisson distrib. int. (see Matlab. code for convg.)
    call random_number(nb_alea)
    DO ii = 1,size(SumPart)-1
        IF (nb_alea .lt. SumPart(1)) THEN
                gen_elec = 0
        ELSE IF ((SumPart(ii).le.nb_alea) .and. (nb_alea .lt. SumPart(ii+1))) THEN
                gen_elec = ii
        END IF
    END DO
END FUNCTION gen_elec


END MODULE poisson_unif
