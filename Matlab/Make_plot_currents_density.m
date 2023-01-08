%% Script to make plots for current density 
% and compare them with and without ion induced emissions

%% load partfiles

ions = espic2dhdf5('H2slanted_self_sustained/Resultats/H2Slanted_1e-12.h5')
ions_less = espic2dhdf5('H2slanted_self_iiee_less/Resultats/H2Slanted_1e-12.h5')

%% Normalised particle number

ions.displayenergy
ions_less.displayenergy

%% Charge over normalised time
%true argument: bool to scale time over coll time

out_less = ions_less.displaycharge(true);
out = ions.displaycharge(true);
step = 5068:length(out.time);
figure
    p(1)= plot(out.time(step), out.charge(step), 'linewidth', 2);
    hold on 
    p(2) = plot(out_less.time, out_less.charge, 'linewidth', 2);
    legend('$q_{iiee}$', '$q_0$' ,'Location','best','Interpreter','latex');
    xlabel('t/\tau_d [-]')
    ylabel('Total Charge [C]')
    set(legend,'FontSize',18);
    set (gca, 'fontsize', 22)

%% surface densities 
figure
%%
ions_less.displaySurfFlux(length(ions_less.t2d))
ions_less.displaySurfFlux(length(ions_less.t2d),1,2)
%%
ions.displaySurfFlux(length(ions.t2d))
ions.displaySurfFlux(length(ions.t2d),1,2)


   
%% Make plots currents evolution
% 1st index is species_id: 1 is for electrons 
% 2nd index is all_cur_id: 1 is to combine electrons and ions

 ions.displaytotcurrevol(length(ions.t2d)-200:5:length(ions.t2d),1,1) 
 
 ions_less.displaytotcurrevol(length(ions_less.t2d)-200:5:length(ions_less.t2d),1,1)

