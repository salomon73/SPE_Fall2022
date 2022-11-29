MODULE polynomial


CONTAINS
REAL FUNCTION eval_polynomial(coefficients, valeur)
    REAL, DIMENSION(:) :: coefficients !< polynomial (e.g fitted yield) coeffs
    REAL :: valeur                     !< point where to evaluate polyn
    INTEGER :: ii

    eval_polynomial = 0
    DO ii=1, size(coefficients)
      eval_polynomial = eval_polynomial+coefficients(ii)*valeur**(ii-1)
    END DO
END FUNCTION eval_polynomial
END MODULE polynomial
