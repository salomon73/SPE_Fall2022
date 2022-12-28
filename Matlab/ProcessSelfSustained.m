%% Processing the data for self sustained iiee in several TREX geometries %%
    addpath(genpath('/home/sguincha/espic2d/matlab'))
    ions = espic2dhdf5('H2Slanted_1e-12.h5') 

%% display current collected in different regions
    ions.displaytotcurrevol(1:10:length(ions.t2d))

%% display surface current density collected in different regions
    ions.displaySurfFlux(length(ions.t2d),1,1)

%% process surf flux for electronic current

    subdiv=1;
    species_id=1; %electrons
    timestep = length(ions.t2d)-1000;
    mflux= ions.Metallicflux(timestep,subdiv);
    lflux= -squeeze(ions.Axialflux(timestep,1,species_id))';
    rflux= squeeze(ions.Axialflux(timestep,length(ions.zgrid),species_id))';

    time=ions.t2d(timestep);

    ids=1:length(mflux);
    P=ions.neutcol.neutdens*ions.kb*300/100;% pressure at room temperature in mbar
        Zf=figure('name','fluxevol');
        linew=5;
        %ions.displaysplbound(gca,1e3);
        contour(ions.zgrid*1e3,ions.rgrid*1e3,ions.geomweight(:,:,1),[0 0],'b-','linewidth',1.5);
        hold on
        for i=1:length(mflux.p)
            x=mflux.p{i}(1,:)*1000;
            y=mflux.p{i}(2,:)*1000;
            y(end)=NaN;
            c=mflux.gamma{i}'*ions.qe/(100^2)/P;
            c(c<=0)=NaN;
            patch(x,y,c,'EdgeColor','interp','LineWidth',linew);
            hold on
        end

    x=ions.zgrid(1)*ones(size(ions.rgrid))*1000;
    y=ions.rgrid*1000;
    y(end)=NaN;
    c=lflux*ions.qe/(100^2)/P;
    c(c<=0)=NaN;
    patch(x,y,c,'EdgeColor','interp','LineWidth',linew);


    x=ions.zgrid(end)*ones(size(ions.rgrid))*1e3;
    y=ions.rgrid*1000;
    y(end)=NaN;
    c=rflux*ions.qe/(100^2)/P;
    c(c<=0)=NaN;
    patch(x,y,c,'EdgeColor','interp','LineWidth',linew);
    
    title(sprintf('t=%4.2f [ns]',time*1e9))

    c=colorbar;
    c.Label.String= 'j\cdotn [A/(cm^2 mbar)]';
    xlabel('z [mm]')
    ylabel('r [mm]')
    colormap(jet)
    set(gca,'colorscale','log')
   
    
    
%% 
    time=ions.t2d(timestep);
    flux     = rflux*ions.qe/(100^2)/P;
    [val Id] = findpeaks(flux,'MinPeakDistance',6);
    peak1 = flux(Id(2));
    peak2 = flux(Id(3));
    
    figure 
    hold on 
    plot(ions.rgrid*1000,flux, 'r-', 'linewidth', 2)
    xlabel('r [mm]', 'interpreter', 'latex')
    ylabel('J', 'interpreter', 'latex')
    set (gca, 'FontSize', 24)
  
%% several points

    for ii=length(ions.t2d)-100:length(ions.t2d)
        mflux= ions.Metallicflux(ii,subdiv);
        lflux= -squeeze(ions.Axialflux(ii,1))';
        rflux= squeeze(ions.Axialflux(ii,length(ions.zgrid)))';

        time=ions.t2d(ii);

        ids=1:length(mflux);
        P=ions.neutcol.neutdens*ions.kb*300/100;% pressure at room temperature in mbar
        flux     = rflux*ions.qe/(100^2)/P;
        [val, Ind] = findpeaks(flux,'NPeaks',2,'SortStr','descend'); %sort(flux);
        peak1(ii-(length(ions.t2d)-100)+1) = flux(Ind(1));
        peak2(ii-(length(ions.t2d)-100)+1) = flux(Ind(2));
%         peak1(ii) = flux(Ind(end));
%         peak2(ii) = flux(Ind(end-1));
    end
%%
    avg1 = smoothdata(peak1);
    figure
    hold on 
    plot(peak1, 'k-')
    hold on 
    plot(peak2, 'r-')
   
 %% several points - all %%
 % good for H2Slanted TREX geometry %
    compteur = 0;
    species_id = 1; % electrons
    for ii=2500:50:length(ions.t2d)
        compteur = compteur +1;
        mflux= ions.Metallicflux(ii,subdiv);
        lflux= -squeeze(ions.Axialflux(ii,1,species_id))';
        rflux= squeeze(ions.Axialflux(ii,length(ions.zgrid), species_id))';

        time=ions.t2d(ii);

        ids=1:length(mflux);
        P=ions.neutcol.neutdens*ions.kb*300/100;% pressure at room temperature in mbar
        flux     = rflux*ions.qe/(100^2)/P;
        [val, Ind] = findpeaks(flux,'NPeaks',2,'SortStr','descend'); %sort(flux);
        peak1(compteur) = flux(Ind(1));
        peak2(compteur) = flux(Ind(2));
        temps(compteur) = time;
%         peak1(ii) = flux(Ind(end));
%         peak2(ii) = flux(Ind(end-1));
    end

%%
    avg1 = smoothdata(peak1);
    avg2 = smoothdata(peak2);
    figure
    hold on 
    p = plot(temps, peak1, 'k-', 'linewidth', 2);
    set(p(1), 'Color', '#D95319')
    hold on 
    p1 = plot(temps, peak2, 'r-', 'linewidth', 2);
    set(p1(1), 'Color', '#0072BD')
    hold on 
    %plot(temps, avg1, 'g--', 'linewidth', 2)
    hold on 
    %plot(temps, avg2, 'b--', 'linewidth', 2)
    xlabel('t [s]', 'interpreter', 'latex')
    ylabel('$J_{max}$ [A/(cm$^2\cdot$mbar)]', 'interpreter', 'latex')
    legend('$J_{max}^{i}$', '$J_{max}^{c}$','Location','best','Interpreter','latex');
    set (gca, 'FontSize', 24)
    set (legend, 'FontSize', 20)
    
    
%% process the surf flux for ionic current

    species_id = 2;
    mflux= ions.MetallicFlux_species(timestep,subdiv);
    lflux= -squeeze(ions.Axialflux(timestep,1,species_id))';
    rflux= squeeze(ions.Axialflux(timestep,length(ions.zgrid),species_id))';
    time=ions.t2d(timestep);
    qe = ions.species(species_id-1).q;
    ids=1:length(mflux);
    
            P=ions.neutcol.neutdens*ions.kb*300/100;% pressure at room temperature in mbar
            f=figure('name','fluxevol');
            linew=5;
            %ions.displaysplbound(gca,1e3);
            contour(ions.zgrid*1e3,ions.rgrid*1e3,ions.geomweight(:,:,1),[0 0],'b-','linewidth',1.5);
            hold on
            for i=2%:length(mflux.p)
                x=mflux.p{i}(1,:)*1000;
                y=mflux.p{i}(2,:)*1000;
                y(end)=NaN;
                %c=mflux.gamma{i}'*ions.qe/(100^2)/P;
                c=mflux.gamma{i}'*qe/(100^2)/P;
                c(c<=0)=NaN;
                patch(x,y,c,'EdgeColor','interp','LineWidth',linew);
                hold on
            end
            
            x=ions.zgrid(1)*ones(size(ions.rgrid))*1000;
            y=ions.rgrid*1000;
            y(end)=NaN;
            %c=lflux*ions.qe/(100^2)/P;
            c=lflux*qe/(100^2)/P;
            c(c<=0)=NaN;
            patch(x,y,c,'EdgeColor','interp','LineWidth',linew);
            
            
            x=ions.zgrid(end)*ones(size(ions.rgrid))*1e3;
            y=ions.rgrid*1000;
            y(end)=NaN;
            %c=rflux*ions.qe/(100^2)/P;
            c=rflux*qe/(100^2)/P;
            c(c<=0)=NaN;
            patch(x,y,c,'EdgeColor','interp','LineWidth',linew);
            
            title(sprintf('t=%4.2f [ns]',time*1e9))
            
            c=colorbar;
            c.Label.String= 'j\cdotn [A/(cm^2 mbar)]';
            xlabel('z [mm]')
            ylabel('r [mm]')
            colormap(jet)
            set(gca,'colorscale','log')