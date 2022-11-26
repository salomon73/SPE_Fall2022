%% Script to find which probability distribution fits better the data %%
%  for electrons energys 

Emin = 0;  % eV
Emax = 25; % eV
x = linspace(Emin,Emax,1000);
theta = 2.0;
kappa = 2.0;

PDF = 1/(gamma(kappa)*theta^kappa) * x.^(kappa-1).*exp(-x/theta);
CDF = 1/(gamma(kappa))*gammainc(x/theta, kappa); 

figure
    subplot(1,2,1)
        hold on 
        plot(x,PDF, 'b-', 'linewidth', 2)
            xlabel('E [eV]', 'Interpreter', 'Latex') 
            ylabel('P(E)', 'Interpreter', 'Latex')
            legend('$\theta = 2.0$ , $\kappa = 2.0$' ,'Location','best','Interpreter','latex');
            set(legend,'FontSize',18);
            set (gca, 'fontsize', 22)
    subplot(1,2,2)
        hold on 
        plot(x,CDF, 'b-', 'linewidth', 2)
            xlabel('E [eV]', 'Interpreter', 'Latex') 
            ylabel('$P(E_{el} \leq E)$', 'Interpreter', 'Latex')
            legend('$\theta = 2.0$ , $\kappa = 2.0$' ,'Location','best','Interpreter','latex');
            set(legend,'FontSize',18);
            set (gca, 'fontsize', 22)

            
%% generates random number follong this gamma distribution %%

nit = 100000;
rand_num = zeros(1,nit);
for ii=1:nit
    beta = rand;
    [val, Ind] = min(abs(beta-CDF));
    rand_num(ii)= x(Ind);
end
mean(rand_num)