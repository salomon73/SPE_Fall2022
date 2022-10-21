
function f=display2Dpotentialwell(obj,timestep,rcoord,clims,rescale, mag,nblv)
            % Display the 2D potential well at time obj.t2d(timestep)
            % if rcoord is true, the potential is evaluated at grid points in r,z coordinates
            % if false, the potential is evaluated at grid points in magnetic field line coordinates
            % clims are the values limits in eV for the color coding of the
            % potential well
            if iscell(timestep)
                timestep=cell2mat(timestep);
            end
            if nargin <3
                rcoord=true;
            end
            if nargin<4
                clims=[-inf inf];
            end
            if nargin<5
                rescale=false;
            end
            if nargin<6 
                mag=[];
            end
            if nargin<7 
                nblv=400;
            end
            
            f=figure('Name',sprintf('%s Potential well',obj.name));
            ax1=gca;
            model=obj.potentialwellmodel(timestep,false,mag,nblv);
            z=model.z;
            r=model.r;
            Pot=model.pot;
            rathet=model.rathet;
            
            N0=obj.N(:,:,1);
            id=find(timestep==0);
            timestep(id)=[];
            Nend=obj.N(:,:,timestep);
            if(~isempty(id))
                N0=zeros(obj.N.nr,obj.N.nz);
                Nend=cat(3,Nend(:,:,1:id-1),N0,Nend(:,:,id:end));
            end
            Nend=mean(Nend,3);
            geomw=obj.geomweight(:,:,1);
            %z(isnan(Pot))=[];
            %r(isnan(Pot))=[];
            %Pot(isnan(Pot))=[];
            if rescale
                if obj.spl_bound.nbsplines >0
                    Pot=Pot/(obj.spl_bound.boundary(2).Dval-obj.spl_bound.boundary(1).Dval);
                else
                    Pot=Pot/(obj.potout-obj.potinn)/obj.phinorm;
                end
            end
            if rcoord
                [Zmesh,Rmesh]=meshgrid(obj.zgrid,obj.rgrid);
                Pot=griddata(z,r,Pot,Zmesh,Rmesh,'natural');
                %Pot(obj.geomweight(:,:,1)<=0)=NaN;
                %Pot=imgaussfilt(Pot,'FilterSize',11); 
                Pot(Pot<=0)=NaN;
                Pot(obj.geomweight(:,:,1)<0)=NaN;
                contourf(obj.zgrid(1:end)*1e3,obj.rgrid(1:end)*1e3,Pot(1:end,1:end),50,'edgecolor','none','Displayname','Well')
                xlabel('z [mm]')
                ylabel('r [mm]')
                xlim([obj.zgrid(1) obj.zgrid(end)]*1e3)
                ylim([obj.rgrid(1) obj.rgrid(end)]*1e3)
                hold(gca, 'on')
                
                rdisp=obj.rgrid;
                
                %% Magnetic field lines
                Blines=obj.rAthet;
                levels=linspace(min(Blines(obj.geomweight(:,:,1)>0)),max(Blines(obj.geomweight(:,:,1)>0)),20);
                Blines(obj.geomweight(:,:,1)<0)=NaN;
                [~,h1]=contour(obj.zgrid*1000,obj.rgrid*1000,Blines,real(levels),'k-.','linewidth',1.2,'Displayname','Magnetic field lines');
                
            else
                rdisp=sort(obj.rAthet(:,end));
                [Zmesh,Rmesh]=meshgrid(obj.zgrid,rdisp);
                Pot=griddata(z,rathet,Pot,Zmesh,Rmesh);
                [Zinit,~]=meshgrid(obj.zgrid,obj.rAthet(:,1));
                %         if  isempty(obj.maxwellsrce)
                %                 end
                N0=griddata(Zinit,obj.rAthet,N0,Zmesh,Rmesh);
                Nend=griddata(Zinit,obj.rAthet,Nend,Zmesh,Rmesh);
                geomw=griddata(Zinit,obj.rAthet,geomw,Zmesh,Rmesh);
                Pot(geomw<=0)=0;
                surface(obj.zgrid(1:end),rdisp,Pot(1:end,1:end),'edgecolor','none')
                ylabel('rA_\theta [Tm^2]')
                xlabel('z [m]')
                xlim([obj.zgrid(1) obj.zgrid(end)])
                ylim([min(rdisp) max(rdisp)])
                hold(gca, 'on')
                
            end
            if (timestep==0)
                title(sprintf('Potential well Vacuum'))
            else
                title(sprintf('Potential well t=%1.2f [ns]',obj.t2d(timestep)*1e9))
            end
            maxdensend=max(Nend(:));
            contourscale=0.1;
            Nend=(Nend-contourscale*maxdensend)/maxdensend;
            maxdens0=max(N0(:));
            contourscale=0.1;
            N0=(N0-contourscale*maxdens0)/maxdens0;
            contour(obj.zgrid*1e3,rdisp*1e3,Nend,linspace(0,1-contourscale,5),'k--','linewidth',1.5,'Displayname','Cloud Boundaries');
            contour(obj.zgrid*1e3,rdisp*1e3,N0,[0 0],'k-.','linewidth',1.5,'Displayname','Source boundaries');
            %contour(obj.zgrid*1e3,rdisp*1e3,geomw,[0 0],'-','linecolor',[0.5 0.5 0.5],'linewidth',1.5,'Displayname','Vessel Boundaries');
            c=colorbar;
            colormap('jet');
            
            
            
            % Grey outline showing the metalic walls
            geomw(obj.geomweight(:,:,1)>0)=-1;
            geomw(obj.geomweight(:,:,1)<=0)=1;
            [c1,hContour]=contourf(ax1,obj.zgrid*1000,obj.rgrid*1000,geomw, [0 0]);
            
            drawnow;
            xlim(ax1,[obj.zgrid(1)*1000 obj.zgrid(end)*1000])
            if(obj.conformgeom)
                ylim([ax1 ],[obj.rgrid(1)*1000 obj.rgrid(rgridend)*1000])
            else
                ylim([ax1],[obj.rgrid(1)*1000 obj.rgrid(end)*1000])
            end
            %ylim(ax1,[0.05*1000 obj.rgrid(end)*1000])
            %xlim([obj.zgrid(1) 0.185]*1e3)
            xlabel(ax1,'z [mm]')
            ylabel(ax1,'r [mm]')
            view(ax1,2)
            if rescale
                c.Label.String= 'well depth / bias';
            else
                c.Label.String= 'well depth [eV]';
            end
            f.PaperUnits='centimeters';
            caxis(clims)
            grid(ax1, 'on');
            hFills=hContour.FacePrims;
            [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
            try
                drawnow
                hFills(1).ColorData = uint8([150;150;150;255]);
                for idx = 2 : numel(hFills)
                    hFills(idx).ColorData(4) = 0;   % default=255
                end
            catch
            end
            
            % add central and external metallic walls if we have a coaxial
            % configuration
            if( obj.walltype >=2 && obj.walltype<=4)
                rectangle('Position',[obj.zgrid(1) obj.r_b obj.zgrid(end)-obj.zgrid(1) 0.001]*1e3,'FaceColor',[150 150 150]/255,'Edgecolor','none')
                ylimits=ylim;
                ylim([ylimits(1),ylimits(2)+1])
            end
            if sum(obj.geomweight(:,1,1))==0
                rectangle('Position',[obj.zgrid(1) obj.r_a-0.001 obj.zgrid(end)-obj.zgrid(1) 0.001]*1e3,'FaceColor',[150 150 150]/255,'Edgecolor','none')
                ylimits=ylim;
                ylim([ylimits(1)-1,ylimits(2)])
            end
            
            %axis equal
            %xlim([-100 200])
            [max_depth,id]=max(abs(Pot(:)));
            [idr,idz]=ind2sub(size(Pot),id);
            fprintf('Maximum potential wel depth: %f eV\n',max_depth)
            fprintf('at location r=%f z=%f [mm]\n',obj.rgrid(idr)*1e3, obj.zgrid(idz)*1e3)
            papsize=[14 8];
            if rcoord
                obj.savegraph(f,sprintf('%s/%s_wellr%i',obj.folder,obj.name,floor(mean(timestep))),papsize);
            else
                obj.savegraph(f,sprintf('%s/%s_wellpsi%i',obj.folder,obj.name,floor(mean(timestep))),papsize);
            end
        end