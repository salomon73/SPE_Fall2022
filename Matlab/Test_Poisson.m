%% -------  Test Poisson distribution -------- %%
%                                               %
% Tests randomly generated int. expected        %
% to follow  a Poisson distr. of parameter      %
% lambda (see function Poisson from Poisson.m)  % 
%                                               %
% S.Guinchard (11/19/2022)                      %
%-----------------------------------------------%

format long               % change double format
it      = 100000;         % iterations #
lambda  = 1;              % Poisson parameter
kmax    = 10;             % highest integ. to be gen. 
vect    = zeros(1,kmax);  % PDF(k) values 
SumPart = zeros(1,kmax);  % CDF(k) values

for ii =  1:length(vect)
   vect(ii) = exp(-lambda)*lambda^(ii-1)/factorial(ii-1); 
   SumPart(ii) = sum(vect(1:ii));
end

% Check that CDF does indeed converge to 1
CumulPoisson = sum(vect, 'all');
disp(strcat('CDF(k=',num2str(kmax),')= ', num2str(CumulPoisson)))

% gen. rand. vals.
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

% # of counts per value k 
[~,~,ix] = unique(I);
C = accumarray(ix,1).'
clear kval
for ii = 1:length(C)    
    kval(ii) = ii-1; % k=0 is an outcome
end

% cat. vector of k to add non-obtained values with count 0
k      = [kval linspace(kval(end)+1, kmax, kmax-kval(end))];
Counts = [C zeros(1,kmax-kval(end))];
SumCumul = [SumPart SumPart(end)+exp(-lambda)*lambda^(kmax)/factorial(kmax)];

% Plot CDF and counts converging to PDF
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


