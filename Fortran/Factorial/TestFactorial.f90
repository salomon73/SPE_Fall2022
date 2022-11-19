PROGRAM test_factorial

USE factorial

INTEGER:: n, resultat 

print*, 'enter an integer n:'
read*,n
resultat = factorial_fun(n)

print*, resultat


END PROGRAM test_factorial 
