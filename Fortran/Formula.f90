! Calculates number of accumulated AIDS cases in USA

PROGRAM FORMULA
        INTEGER T
        REAL VAL
        READ*, T
        VAL = 174.6 * (T - 1981.2) ** 3
        PRINT*, 'Accumulated cases in by unit time', T, ':', VAL
END PROGRAM FORMULA

