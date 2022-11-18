%% Test Poisson distribution %%
% from generation of a random number between min and max 
% of cumulative random distribution 

% Redo the code correctly and show that it converges indeed towards a
% poisson distribution
format long
it      = 1000000;
lambda  = 0.25;
kmax    = 15;
vect    = zeros(1,kmax);
SumPart = zeros(1,kmax); 

for ii =  1:length(vect)
   vect(ii) = exp(-lambda)*lambda^(ii-1)/factorial(ii-1); 
end
for jj = 1:length(SumPart)
    SumPart(jj) = sum(vect(1:jj));
end
CumulPoisson = sum(vect, 'all');


%% 
for ii =1: it 
   alea   = rand;
   for jj=1:length(SumPart)-1 
       if(alea< SumPart(1))
            I(ii)=0;
       elseif(le(SumPart(jj), alea) && alea < SumPart(jj+1))
            I(ii) = jj;
       end
   end
    
end


%%
[~,~,ix] = unique(I);
C = accumarray(ix,1).'
clear kval
for ii = 1:length(C)    
    kval(ii) = ii-1;
end



k      = [kval linspace(kmax-kval(end), kmax, kmax-kval(end))];
Counts = [C zeros(1,kmax-kval(end))];

SumCumul = [SumPart SumPart(end)+exp(-lambda)*lambda^(15)/factorial(15)];

figure
hold on 
plot(k,1/it*Counts, 'r-', 'linewidth',2)
xlabel('k', 'Interpreter', 'Latex')
ylabel('counts', 'Interpreter', 'Latex')
set (gca, 'fontsize', 20) 

figure
plot(k,SumCumul, 'r-', 'linewidth',2)
xlabel('k', 'Interpreter', 'Latex')
ylabel('CDF', 'Interpreter', 'Latex')
set (gca, 'fontsize', 20) 


