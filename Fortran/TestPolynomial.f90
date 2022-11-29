PROGRAM test_polynomial

        USE polynomial

        REAL, DIMENSION(2) :: coefficients_1H  = (/0.0, 2.9778E02 /)
        REAL, DIMENSION(2) :: coefficients_1He = (/0.1273010048, 1.70478995200000E02 /)
        REAL, DIMENSION(2) :: coefficients_1Ne = (/0.0518524160, 2.45927583999999E02 /)
        real ::energy = 0
        integer :: neuttype_id = 0
        real :: eval_pol_val = 0
        print*, "enter energy value"
        read*, energy
        print*, "select neuttype"
        read*, neuttype_id
        select case(neuttype_id)
                case(1)
                        eval_pol_val = eval_polynomial(coefficients_1H, energy)

                case(2)
                        eval_pol_val = eval_polynomial(coefficients_1He, energy)

                case(3)
                        eval_pol_val = eval_polynomial(coefficients_1Ne, energy)
                case default 
                        eval_pol_val = eval_polynomial(coefficients_1H, energy)

        end select         
print*, eval_pol_val

END PROGRAM test_polynomial
