!------------------------------------------------------------------------------
! EPFL/Swiss Plasma Center
!------------------------------------------------------------------------------
!
! MODULE: materials 
!
!> @author / date
!> S. Guinchard - EPFL/SPC / 11-25-22
!
! DESCRIPTION:
!> Module storing coefficients corresponding to polynomial
!> fits of energy loss function dE/dx in several materials
!> N.B. These coefficients fit dE/dx in MeV/cm for protons 
!> only and should be adapted for other ion species
!------------------------------------------------------------------------------

MODULE materials 

        USE constants 

        IMPLICIT NONE    
     ! ----------------- 304 STAINLESS STEEL ----------------
     !> All the coefficients below were obtained by piecewise
     !> fit of dE/dx curve for 304 stainless steel
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1H_SS  = (/0.0, 2.9778E02 /)
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1He_SS = (/0.1273010048, 1.70478995200000E02 /)
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1Ne_SS = (/0.0518524160, 2.45927583999999E02 /)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_2_SS = &
        (/1.834520315818557E02,1.320304216355084E05,-8.583700602370013E06,3.526140145560557E08/)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_3_SS = &
        (/2.471679999999999E02,9.695466666666670E04,-2.475200000000003E06, 2.325333333333340E07/)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_4_SS = &
        (/2.533904454349683E03,-1.766016382825937E05,8.202640024592019E06, -1.125320217235288E08/)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_5_SS = &
        (/8.786057142856745E02,2.856595238095524e+04,-1.834285714286391E05, 3.333333333338620E05/)

     ! --------------------- COPPER -------------------------
     !> All the coefficients below were obtained by piecewise
     !> fit of dE/dx curve for Copper 
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1H_Cu  = (/0.0, 2.812300E02 /)
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1He_Cu = (/0.1273010048, 1.539289952000001E02 /)
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1Ne_Cu = (/0.0518524160, 2.293775840000000E02 /)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_2_Cu = &
        (/1.703970762187365E02, 1.273920899063999E05, -8.680492697427949E06, 3.535110244304671E08/)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_3_Cu = &
        (/2.618095714285664E02, 8.520814285714373E04, -2.094971428571471E06, 2.064000000000054E07/)   
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_4_Cu = &
        (/1.030650000000107E03, -6.231904761918026E03, 1.468571428571969E06, -2.506666666667393E07/)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_5_Cu = &
        (/7.117599999999986E02, 3.874166666666669E04, -5.099999999999992E05, 2.733333333333320E06/)

     ! --------------------- ALUMINUM ----------------------
     !> All the coefficients below were obtained by piecewise
     !> fit of dE/dx curve for Aluminum
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1H_Al  = (/0.0, 2.100900000000000E02 /)
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1He_Al = (/0.1273010048, 82.788995200000016 /)
        REAL(KIND = db), DIMENSION(2), PARAMETER  :: coefficients_1Ne_Al = (/0.0518524160, 1.582375839999999E02 /)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_2_Al = &
        (/1.224483183711612E02, 9.731277035468460E04, -6.878546755544925E06, 2.795627748449915E08 /)    
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_3_Al = &
        (/1.998525714285755E02, 6.229714285714197E04, -1.531771428571365E06, 1.567999999999851E07/)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_4_Al = &
        (/8.135440000000424E02, -1.333600000000502E04, 1.553600000000198E06, -2.624000000000260E07/)
        REAL(KIND = db), DIMENSION(4), PARAMETER  :: coefficients_5_Al = &
        (/4.052485714285687E02, 3.974642857142875E04, -6.631428571428614E05, 3.800000000000037E06 /)



END MODULE materials
