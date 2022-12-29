%% physical parameters
timestep = length(ions.t2d)-100:10:length(ions.t2d);
mfluxe = ions.Metallicflux(timestep,subdiv);
mfluxi = ions.MetallicFlux_species(timestep,subdiv);
qe = abs(ions.qe);
P=ions.neutcol.neutdens*obj.kb*300/100;
%% plot electronic current flux
figure
hold on 
plot(mfluxe.p{1}(1,:), qe/P*mfluxe.gamma{1}(:,1)')
xlabel('Z [m]', 'interpreter', 'latex')
ylabel('$\Phi_e$ [A/m$^2$/mbar]', 'interpreter', 'latex')
set (gca, 'fontsize', 20)

%% plot ionic current flux
figure
hold on 
plot(mfluxi.p{2}(1,:), qe/P*mfluxi.gamma{2}(:,1)')
xlabel('Z [m]', 'interpreter', 'latex')
ylabel('$\Phi_i$ [A/m$^2$/mbar]', 'interpreter', 'latex')
set (gca, 'fontsize', 20)

%% Currents 
qe = abs(ions.qe);
for ii =1:length(mfluxi.gamma) 
    for jj = 1:length(timestep)

        flux = (qe/P)*mfluxi.gamma{ii}(:,jj)'.*mfluxi.p{ii}(2,:);
        flux = flux(~isnan(flux));
        current(ii,jj) =  2*pi*trapz(mfluxi.p{ii}(1,~isnan(flux)), flux);  
    end
end
%%
qe = ions.qsim/ions.weight;
for ii = 1:length(mfluxe.gamma) 
    for jj = 1:length(timestep)

        flux = (qe/P)*mfluxe.gamma{ii}(:,jj)'.*mfluxe.p{ii}(2,:);
        flux = flux(~isnan(flux));
        current(ii,jj) =  2*pi*trapz(mfluxe.p{ii}(1,~isnan(flux)), flux);  
    end
end
%%
figure
hold on 
plot(current(1,:))
hold on 
plot(current(2,:))