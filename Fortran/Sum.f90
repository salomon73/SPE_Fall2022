PROGRAM sommation
        
        INTEGER, DIMENSION(3):: tab, diff
        REAL :: resultat
        INTEGER, DIMENSION(10) :: vecteur 
       
       
       
        tab = (/1,2,3/)
        resultat = sum(tab)
        print*, resultat
        resultat = sum(tab(1:2))
        print*, resultat

        DO i=1,10
          vecteur(i) = 1 
        END DO 
        PRINT*, vecteur
        diff = abs(tab-tab)

        PRINT*, diff



END PROGRAM sommation



