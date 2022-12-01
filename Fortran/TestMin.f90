PROGRAM test_min_array

        REAL, DIMENSION(3) :: tableau, tableau2, tableau3, Diff
        INTEGER :: ii, posid
        REAL :: valeur
        !DO ii=1,size(tableau)
        !        tableau(ii) = 4-ii
        !END DO 
        tableau = (/2.0,1.0,4.0/)
        tableau2 = 3* tableau
        !do ii = 1,size(tableau)
        Diff = 3-tableau
        !end do
        print*, tableau
        print*, tableau2
        print*, Diff 
        posid = minloc(tableau,dim = 1)
        valeur   = minval(tableau, dim=1)
        print*, posid
        print*, valeur

END PROGRAM test_min_array
