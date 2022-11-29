%% Script to find which probability distribution fits better the data %%
%  for electrons energys 

Emin = 0;  % eV
Emax = 25; % eV
x = linspace(Emin,Emax,1000);
theta = .5;
kappa = 4.0;

PDF1 = 1/(gamma(kappa)*theta^kappa) * x.^(kappa-1).*exp(-x/theta);
CDF1 = gammainc(x/theta,kappa); 

theta2 = 2.0;
kappa2 = 1.0;


PDF2 = 1/(gamma(kappa2)*theta2^kappa2) * x.^(kappa2-1).*exp(-x/theta2);
CDF2 = gammainc(x/theta2,kappa2); 

theta3 = 2/3;
kappa3 = 3.0;


PDF3 = 1/(gamma(kappa3)*theta3^kappa3) * x.^(kappa3-1).*exp(-x/theta3);
CDF3 = gammainc(x/theta3,kappa3); 

figure
    subplot(1,2,1)
        hold on 
        plot(x,PDF1, 'b-', 'linewidth', 2)
        hold on 
        plot(x,PDF2, 'r-', 'linewidth', 2)
        hold on 
        plot(x,PDF3, 'k-', 'linewidth', 2)
            xlabel('E [eV]', 'Interpreter', 'Latex') 
            ylabel('P(E)', 'Interpreter', 'Latex')
            legend(strcat('$\theta =$', num2str(theta),'$\kappa =$', num2str(kappa)) ,...
                   strcat('$\theta =$', num2str(theta2),'$\kappa =$', num2str(kappa2)),...
                   strcat('$\theta =$', num2str(theta3),'$\kappa =$', num2str(kappa3))...
                   ,'Location','best','Interpreter','latex');
            set(legend,'FontSize',18);
            set (gca, 'fontsize', 22)
    subplot(1,2,2)
        hold on 
        plot(x,CDF1, 'b-', 'linewidth', 2)
        hold on 
        plot(x,CDF2, 'r-', 'linewidth', 2)
        hold on 
        plot(x,CDF3, 'k-', 'linewidth', 2)
            xlabel('E [eV]', 'Interpreter', 'Latex') 
            ylabel('$P(E_{el} \leq E)$', 'Interpreter', 'Latex')
            legend(strcat('$\theta =$', num2str(theta),'$\kappa =$', num2str(kappa)) ,...
                   strcat('$\theta =$', num2str(theta2),'$\kappa =$', num2str(kappa2)),...
                   strcat('$\theta =$', num2str(theta3),'$\kappa =$', num2str(kappa3))...
                   ,'Location','best','Interpreter','latex');
            set(legend,'FontSize',18);
            set (gca, 'fontsize', 22)

      %strcat('$\theta =$', num2str(theta), {'  '}  ,'$\kappa =$', num2str(kappa))      
%% generates random number follong this gamma distribution %%

nit = 100000;
rand_num = zeros(1,nit);
for ii=1:nit
    beta = rand;
    [val, Ind] = min(abs(beta-CDF1));
    rand_num(ii)= x(Ind);
%     num = gamrnd(2.0,2.0);
%     rand_num(ii) = num;
end
mean(rand_num)