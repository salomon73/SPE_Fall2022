function displaySurfFluxes(obj,timesteps, ids)
            %displaySurfFluxes plot the time evolution of the current
            %densities on the domain boundaries for times t2d(timesteps)
            %also plot the different boundaries to see which flux belong to
            %which boundary
            
            mflux= obj.Metallicflux(timesteps);
            lflux= -squeeze(obj.Axialflux(timesteps,1))';
            rflux= squeeze(obj.Axialflux(timesteps,length(obj.zgrid)))';
            
            time=obj.t2d(timesteps);
            
            if nargin<3
                ids=1:(size(mflux.p,2)+2);
            end
            
            %%
            P=obj.neutcol.neutdens*obj.kb*300/100;% pressure at room temperature in mbar
            f=figure('name','fluxevol');
            tiledlayout('flow')
            j=1;
            for i=1:length(mflux.p)
                if(find(i+2==ids))
                    ax(j)=nexttile;
                    j=j+1;
                    if  issorted(mflux.p{i}(1,:),'strictascend')
                        contourf(mflux.p{i}(1,:)*100,time*1e9,mflux.gamma{i}'*obj.qe/(100^2)/P,'linestyle','none')
                        xlabel('z [cm]')
                    else
                        contourf(linspace(0,1,length(mflux.p{i}(1,:))),time*1e9,mflux.gamma{i}'*obj.qe/(100^2)/P,'linestyle','none')
                        xlabel('s [-]')
                    end
                    title(sprintf('Wall %i',i))
                    
                    c=colorbar;
                    c.Label.String= 'j\cdotn [A/(cm^2 mbar)]';
                end
            end
            
            if(find(1==ids))
                ax(j)=nexttile;
                contourf(obj.rgrid*100,time*1e9,lflux*obj.qe/(100^2)/P,'linestyle','none')
                title('left')
                xlabel('r [cm]')
                c=colorbar;
                c.Label.String= 'j\cdotn [A/(cm^2 mbar)]';
            end
            
            if(find(2==ids))
                ax(j+1)=nexttile;
                contourf(obj.rgrid*100,time*1e9,rflux*obj.qe/(100^2)/P,'linestyle','none')
                title('right')
                xlabel('r [cm]')
                c=colorbar;
                c.Label.String= 'j\cdotn [A/(cm^2 mbar)]';
            end
            
            
            ylabel(ax,'t [ns]')
            
            nexttile;
            for i=1:length(mflux.p)
                plot(mflux.p{i}(1,:)*100,mflux.p{i}(2,:)*100,'displayname',sprintf('Wall %i',i),'linewidth',2)
                hold on
            end
            legend('location','eastoutside')
            title('Domain')
            plot(ones(size(obj.rgrid))*obj.zgrid(1)*100,obj.rgrid*100,'displayname','left','linewidth',2)
            plot(ones(size(obj.rgrid))*obj.zgrid(end)*100,obj.rgrid*100,'displayname','right','linewidth',2)
            xlabel('z [cm]')
            ylabel('r [cm]')
            %xlim([obj.zgrid(1) obj.zgrid(end)]*100)
            %ylim([obj.rgrid(1) obj.rgrid(end)]*100)
            
            obj.savegraph(f,sprintf('%s/%s_surfFluxEvol',obj.folder,obj.name),[16 14]);
        end