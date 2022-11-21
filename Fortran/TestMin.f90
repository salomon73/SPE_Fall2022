PROGRAM test_min_array

        REAL, DIMENSION(3) :: tableau
        INTEGER :: ii
        REAL :: posid, valeur
        !DO ii=1,size(tableau)
        !        tableau(ii) = 4-ii
        !END DO 
        tableau = (/2.0,1.0,4.0/)

        print*, tableau
        posid = minloc(tableau,dim = 1)
        valeur   = minval(tableau, dim=1)
        print*, posid
        print*, valeur

END PROGRAM test_min_array
