%% Script to find which probability distribution fits better the data %%
%  for electrons energys 

Emin = 0;  % eV
Emax = 25; % eV
x = linspace(Emin,Emax,500);
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
        h1 = plot(x,PDF1, '-', 'linewidth', 2);
        set(h1(1), 'Color', '#0072BD')
        hold on 
        h2 = plot(x,PDF2, '-'	, 'linewidth', 2);
        set(h2(1), 'Color', '#D95319')
        hold on 
        h3 = plot(x,PDF3 ,'-', 'linewidth', 2);
        set(h3(1), 'Color', '#77AC30')
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
        h1 = plot(x,CDF1, '-', 'linewidth', 2);
        set(h1(1), 'Color', '#0072BD')
        hold on 
        h2 = plot(x,CDF2, '-'	, 'linewidth', 2);
        set(h2(1), 'Color', '#D95319')
        hold on 
        h3 = plot(x,CDF3 ,'-', 'linewidth', 2);
        set(h3(1), 'Color', '#77AC30')
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
end
mean(rand_num)



N = hist(rand_num)
figure
    histogram(rand_num, 'Normalization', 'PDF', 'FaceColor', '#0072BD', 'FaceAlpha', 0.5)
    hold on 
    h1 = plot(x,PDF1, 'b-', 'linewidth', 2)
    set(h1(1), 'Color', '#D95319')
    xlabel('E [eV]', 'Interpreter', 'Latex') 
    ylabel('$P(E_{el} = E)$', 'Interpreter', 'Latex')
    legend( '$\frac{n_{counts}}{n_{it}}(k)$','PDF(E)', ...
       'Location','northeast','Interpreter','latex');
    set(legend,'FontSize',20);
    set (gca, 'fontsize', 22)
