%%
    timesteps=1:length(ions.t2d);
    species_id =1;
    toptitle="";
    scalet=true;
    dens=true;
    subdiv=1;
    nmean=1;
%%
        if scalet
            if ions.neutcol.present
                vexb0=(ions.Ez(:,:,1).*ions.Br'-ions.Er(:,:,1).*ions.Bz')./(ions.B'.^2);
                vexb0(ions.geomweight(:,:,1)<=0)=0;
                E=0.5*ions.msim/ions.weight*mean(abs(vexb0(:)))^2/ions.qe;
                taucol=1/(ions.neutcol.neutdens*mean(abs(vexb0(:)))*(ions.sigio(E)+ions.sigmela(E)+ions.sigmio(E)));
                try
                    Sio_S=1e17*(ions.neutcol.neutdens*mean(abs(vexb0(:)))*ions.sigio(E))/(ions.maxwellsrce.frequency*ions.weight...
                          /(pi*(ions.maxwellsrce.rlim(2)^2-ions.maxwellsrce.rlim(1)^2)*diff(ions.maxwellsrce.zlim)))
                catch
                end
                tlabel='t/\tau_d [-]';
            else
                taucol=2*pi/ions.omece;
                tlabel='t/\tau_ce [-]';
            end
        else
            taucol=1e-9;
            tlabel='t [ns]';
        end
%%
        if dens
            if species_id==1
                N=ions.N(:,:,timesteps);
            else 
                N=ions.species_moments.N(:,:,timesteps);
            end
            geomw=ions.geomweight(:,:,1);
                geomw(geomw<0)=0;
                geomw(geomw>0)=1;
                N=N.*geomw;
                %[~,idl]=max(N(:,:,end),[],'all','linear');
                %[ir,iz]=ind2sub(size(geomw),idl);
                %nmax=squeeze(max(max(N,[],1),[],2));
                tn=(ions.t2d(timesteps));
                nmax=zeros(2,length(tn));
                nrhalf= floor(0.5*length(ions.rgrid));%find(ions.rgrid>0.07);

                nmax(1,:)=squeeze(max(max(N(1:nrhalf,:,:),[],1),[],2));
                nmax(2,:)=squeeze(max(max(N(nrhalf+1:end,:,:),[],1),[],2));


                nlabel='n_{e,max} [m^{-3}]';
                ndlabel='n_{e,max}';
        else
            t0dst=find(ions.t0d>=ions.t2d(timesteps(1)),1,'first');
            t0dlst=find(ions.t0d<=ions.t2d(timesteps(end)),1,'last');
            tn=ions.t0d(t0dst:t0dlst);
            nmax=ions.npart(t0dst:t0dlst)*ions.weight;
            nlabel='Nb e^-';
            ndlabel='Nb e^-';
        end

        %%
        if species_id ==1
        [currents,pos]=ions.OutCurrents(timesteps,subdiv);
        P=ions.neutcol.neutdens*ions.kb*300/100;% pressure at room temperature in mbar
        currents=currents/P;
        else % for ionic currents
        [currents,pos]=ions.OutCurrents_species(timesteps,subdiv);
        P=ions.neutcol.neutdens*ions.kb*300/100;% pressure at room temperature in mbar
        currents=currents/P;   
        end
        f=figure('Name',sprintf('%s Charges',ions.name));
        tiledlayout(2,1)
        nexttile
        % Plot the evolution of nb of particles
        yyaxis right
        for i=1:size(nmax,1)
            p=plot(tn/taucol,nmax(i,:),'b','linewidth',2.2,'Displayname',sprintf('%s, %d',ndlabel,i));
            hold on
        end
        ylabel(nlabel)
        axl=gca;
        axl.YAxis(2).Color=p.Color;
        ylim([0 inf])

        if(ions.B(1,1)>ions.B(end,1))
            lname='HFS';
            rname='LFS';
        else
            lname='LFS';
            rname='HFS';
        end

        yyaxis( 'left');
        map=colormap(lines);
        set(axl,'linestyleorder',{'-',':','--','*','+'},...
            'ColorOrder',map(2:7,:), 'NextPlot','replacechildren')
        p(1)=plot(axl,ions.t2d(timesteps)/taucol,movmean(currents(1,:),nmean),'Displayname',lname,'linewidth',1.8);
        hold on
        p(2)=plot(axl,ions.t2d(timesteps)/taucol,movmean(currents(2,:),nmean),'Displayname',rname,'linewidth',1.8);
        % Plot the currents
        for i=3:size(currents,1)
            p(i)=plot(axl,ions.t2d(timesteps)/taucol,movmean(currents(i,:),nmean),'Displayname',sprintf('border %i',i-2),'linewidth',1.8);
        end
        plot(axl,ions.t2d(timesteps)/taucol,movmean(sum(currents(:,:),1,'omitnan'),nmean),'k-','Displayname','total','linewidth',1.8);
        hold on 
        xlabel(tlabel)
        ylabel('I/p_n [A/mbar]')
        grid on
        set(gca,'fontsize',12)
        ax.YAxis(1).Color='black';

        %legend('Orientation','horizontal','location','north','numcolumns',3)

        if ~isempty(toptitle)
            title(toptitle)
        end

        ax2=nexttile;
        geomw=ions.geomweight(:,:,1);
        geomw(geomw<=0)=0;
        geomw(geomw>0)=NaN;
        [c1,hContour]=contourf(ax2,ions.zgrid*1000,ions.rgrid*1000,geomw, [0 0]);
        hold on
        drawnow;
        grid on;

        for i=1:length(pos)
            plot(ax2,pos{i}(1,:)*1000,pos{i}(2,:)*1000,'linestyle',p(i+2).LineStyle,...
                'color',p(i+2).Color,'marker',p(i+2).Marker,...
                'displayname',sprintf('border %i',i),'linewidth',1.8)
            hold on
        end
        title('Domain')
        plot(ax2,ones(size(ions.rgrid))*ions.zgrid(1)*1000,ions.rgrid*1000,'linestyle',p(1).LineStyle,...
            'color',p(1).Color,'marker',p(1).Marker,...
            'displayname',lname,'linewidth',1.8)
        plot(ax2,ones(size(ions.rgrid))*ions.zgrid(end)*1000,ions.rgrid*1000,'linestyle',p(2).LineStyle,...
            'color',p(2).Color,'marker',p(2).Marker,...
            'displayname',rname,'linewidth',1.8)
        xlabel('z [mm]')
        ylabel('r [mm]')
        grid on
        set(gca,'fontsize',12)
        hFills=hContour.FacePrims;
        [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
        try
            hFills(1).ColorData = uint8([150;150;150;255]);
            for idx = 2 : numel(hFills)
                hFills(idx).ColorData(4) = 0;   % default=255
            end
        catch
        end
        %legend('Orientation','horizontal','location','north','numcolumns',4)


        fprintf('mean total current: %f [A/mbar]\n',mean(sum(currents(:,max(1,size(currents,2)-30):end),1,'omitnan')));