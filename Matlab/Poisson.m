%% Test Poisson distribution %%
% from generation of a random number between min and max 
% of cumulative random distribution 

% Redo the code correctly and show that it converges indeed towards a
% poisson distribution

it      = 1000000;
lambda  = 1;
kmax    = 10;
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
for ii = 1:it
    
    alea   = rand*ones(1,kmax);
    diff   = abs(alea-SumPart);
    [~,Id] = min(diff);
    I(ii) = Id-1;
end





%%
[~,~,ix] = unique(I);
C = accumarray(ix,1).'
clear kval
for ii = 1:length(C)    
    kval(ii) = ii-1;
end
if(lambda == 1)
    C(2) = 2*C(2);
end
figure
hold on 
plot(kval,C, 'k-', 'linewidth',2)
xlabel('k')
ylabel('counts')
set (gca, 'fontsize', 20) 

figure
plot(SumPart, 'k-', 'linewidth',2)