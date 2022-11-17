PROGRAM POISSON

        USE random_distr
        implicit none

        REAL :: lambda = 0.7
        LOGICAL :: first = .false.
        INTEGER :: num 

        num = 1
        num = random_poisson(lambda, first)
        print*, num


END PROGRAM POISSON 
