MODULE factorial 
CONTAINS

INTEGER FUNCTION factorial_fun(n)
    INTEGER :: n 
    factorial_fun = 1
    IF (n .lt. 0) THEN 
        RETURN 
    ELSE IF (n == 0) THEN 
        factorial_fun = 1
    ELSE 
        DO ii=1,n
         factorial_fun = factorial_fun*ii
        END DO  
   END IF
END FUNCTION factorial_fun


END MODULE factorial 

