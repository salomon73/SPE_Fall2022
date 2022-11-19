PROGRAM POISSON


        USE poisson_unif
        implicit none

        REAL :: lambda = 1
        INTEGER :: num, kmax
        kmax = 15

        num = gen_elec(lambda, kmax)
        print*, num 

END PROGRAM POISSON 
