function k = Poisson(lambda, kmax)
 
 % POISSON( LAMBDA, KMAX)
 %==================
 %
 % Function generating Poisson distributed random integers
 % for a given parameter lambda through randomly generated
 % numbers uniformly between [0,1]
 %
 %
 % INPUT
 % -----
 %   -lambda : Poisson parameter
 %   - kmax  : highest int. to be gen.
 %
 % Written by S.Guinchard (11/19/2022)
 % -----------------------------------
 

 vect    = zeros(1,kmax);  % PDF(k) values 
 SumPart = zeros(1,kmax);  % CDF(k) values
 for ii =  1:length(vect)
   vect(ii) = exp(-lambda)*lambda^(ii-1)/factorial(ii-1); 
   SumPart(ii) = sum(vect(1:ii));
 end
 % Check that CDF does indeed converge to 1
 CumulPoisson = sum(vect, 'all');
 disp(strcat('CDF(k=',num2str(kmax),')= ', num2str(CumulPoisson)))
 if CumulPoisson < 0.994
   warning('CDF(kmax) < 1 : kmax might be too low') 
 end

 alea   = rand;
 for jj=1:length(SumPart)-1 
     if(alea< SumPart(1))
          k=0;
     elseif(le(SumPart(jj), alea) && alea < SumPart(jj+1))
            k = jj;
     end
 end

end