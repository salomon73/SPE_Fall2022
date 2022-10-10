nbzid=15;
zindices=linspace(length(M.zgrid)/6,length(M.zgrid)*5/6,nbzid);
zconsid=1:nbzid;

for zid=zconsid%1:length(zindices)%;%floor(length(M.zgrid)/2);
    zindex=floor(zindices(zid));
    zfolder=sprintf('diocz_%03d',zindex);
    if ~isfolder(zfolder) || ~isfile(sprintf('%s/results.mat',zfolder))
        fprintf('Warning: result %s doesn''t exist\n',zfolder)
        continue
    end
    result=load(sprintf('%s/results.mat',zfolder));
    if ~isfield(result,'geomweight')
        result.geomweight=M.geomweight(:,:,1);
    end
    rindex=M.geomweight(:,zindex,1)>=0;
    
    for l=result.lsteps
        
        r=result.r(:,l);
        phin=result.phinum(:,:,l);
        w=result.w(:,l);
        n=result.n;
        omegar=result.omegadrift;
        rgrid=result.rgrid;
        omegar=omegar(rindex,zindex);
        dn=(n(3:end)-n(1:end-2))./(rgrid(3:end)-rgrid(1:end-2));
        a=min(rgrid);
        b=max(rgrid);
        idmin=find(r==a);
        idmax=find(r==b);
        wd=(M.qsim/M.weight)^2*n/M.eps_0/(M.msim/M.weight)/2/result.wce;
        f=figure;
        subplot(4,1,1)
        hold on
        xlabel('R [m]')
        ylabel('n [m^{-3}]')
        plots=gobjects(2,1);
        plots(1)=plot(rgrid,n,'displayname','n');
        yyaxis right
        plots(2)=plot(rgrid,omegar,'displayname','\omega_{re}');
        hold on
        %plots(3)=plot(rgrid,wd,'displayname','\omega_{d}');
        legend(plots,{'n','\omega_r','\omega_d'},'location','northwest')
        ylabel('\omega_{re} [rad\cdots^{-1}]')
        grid minor
        
        ax=gca;
        ax.LineStyleOrder={'-','--','-.',':'};
        subplot(4,1,[2,3])
        hold on
        p=gobjects();
        j=1;
        for i=1:size(phin,1)
            if imag(w(i))/abs(real(w(i)))>1e-4
                freq=w(i);
                disphinum=phin(i,:);
                disphinum=disphinum/max(abs(disphinum));
                yyaxis left
                p(j)=plot(r(idmin:idmax),[0 real(disphinum(idmin:idmax-2)) 0],'displayname',sprintf('\\omega=%1.3g%+1.3gi',real(freq),imag(freq)));
                j=j+1;
                yyaxis right
                plot(r(idmin:idmax),[0 imag(disphinum(idmin:idmax-2)) 0],'displayname',sprintf('\\omega=%1.3g%+1.3gi',real(freq),imag(freq)))
            end
        end
        subplot(4,1,1)
        if (abs(imag(w(1))/real(w(1)))>1e-10)
            title(sprintf('Unstable case l=%d, \\omega_d=%1.3g, z=%1.3f [mm]',l,max(wd),M.zgrid(zindex)*1e3))
        else
            title(sprintf('Stable case l=%d, \\omega_d=%1.3g z=%1.3f [mm]',l,max(wd),M.zgrid(zindex)*1e3))
        end
        subplot(4,1,[2,3])
        yyaxis left
        ylim([-1 1])
        xlabel('R [m]')
        ylabel('real(\delta\phi) [a.u.]')
        grid on
        yyaxis right
        ylim([-1 1])
        grid minor
        ylabel('imag(\delta\phi) [a.u.]')
        f.PaperUnits='centimeters';
        f.PaperSize=[16,20];
        legend('location','northwest')
        if length(p)>=1 && all(ishghandle(p))
            legend(p)
        end
        subplot(4,1,4)
        plot(real(w),imag(w),'x')
        xlabel('\omega_r [rad\cdots^{-1}]')
        ylabel('\omega_i [rad\cdots^{-1}]')
        subplot(4,1,1)
        ax1=gca;
        subplot(4,1,[2,3])
        ax2=gca;
        drawnow
        posit1=get(ax1,'position');
        posit2=get(ax2,'position');
        posit2(3)=posit1(3);
        set(ax2,'position',posit2);
        print(f,sprintf('%s/%s_dioc_l%d',zfolder,M.name,l),'-dpdf','-fillpage')
        savefig(f,sprintf('%s/%s_dioc_l%d',zfolder,M.name,l))
        close(f) 
        f=figure;
        subplot(2,1,1)
        hold on
        plot(rgrid,n/max(n),'displayname','n')
        plot(rgrid(2:end-1),dn/max(abs(dn)),'displayname','dn')
        plot(rgrid,(real(w(1))-l*omegar)/max(abs((real(w(1))-l*omegar))),'displayname','w-lwre')
        legend
        grid on
        xlabel('r [m]')
        ylabel('Normalized quantity a.u.')
        subplot(2,1,2)
        disphinum=phin(1,:);
        disphinum=disphinum/max(abs(disphinum));
        p(j)=plot(r(idmin:idmax),[0 real(disphinum(idmin:idmax-2)) 0],'displayname',sprintf('\\omega=%1.3g%+1.3gi',real(freq),imag(freq)));
        xlabel('r [m]')
        ylabel('Re(\delta\phi) a.u.')
        f.PaperUnits='centimeters';
        f.PaperSize=[12,6];
        print(f,sprintf('%s/%s_maxw_l%d',zfolder,M.name,l),'-dpdf','-fillpage')
        savefig(f,sprintf('%s/%s_maxw_l%d',zfolder,M.name,l))
        close(f) 
    end
end