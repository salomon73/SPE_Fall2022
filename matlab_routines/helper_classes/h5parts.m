classdef h5parts
    properties
        fullpath
        R
        Z
        THET
        VR
        VZ
        VTHET
        Rindex
        Zindex
        partindex
        partepot
        nlclassical
        q
        m
        weight
        tpart
        it2
        nbparts
        
        rgrid
        zgrid
        vnorm
        rindex
        zindex
        
        file
        
        parent
    
    end
    methods
        function obj=h5parts(filename,hdf5group,parent)
            obj.fullpath=filename;
            obj.file=filename;
            obj.nlclassical=parent.nlclassical;
            obj.it2=parent.it2;
            obj.rgrid=parent.rgrid;
            obj.zgrid=parent.zgrid;
            obj.vnorm=parent.vnorm;
            obj.tpart = h5read(obj.fullpath,sprintf('%s/time',hdf5group));
            obj.nbparts = h5read(obj.fullpath,sprintf('%s/Nparts',hdf5group));
            obj.parent=parent;
            if sum(obj.nbparts)>0
            
            obj.R = h5partsquantity(obj.fullpath,hdf5group,'R');
            obj.Z = h5partsquantity(obj.fullpath,hdf5group,'Z');
            obj.THET = h5partsquantity(obj.fullpath,hdf5group,'THET');
            obj.VR = h5partsquantity(obj.fullpath,hdf5group,'UR',obj.vnorm);
            obj.VZ = h5partsquantity(obj.fullpath,hdf5group,'UZ',obj.vnorm);
            obj.VTHET= h5partsquantity(obj.fullpath,hdf5group,'UTHET',obj.vnorm);
            obj.q= h5readatt(obj.fullpath,hdf5group,'q');
            obj.m= h5readatt(obj.fullpath,hdf5group,'m');
            obj.weight= h5readatt(obj.fullpath,hdf5group,'weight');
            obj.partepot = h5partsquantity(filename,hdf5group,'pot',obj.q);
            obj.rindex=1:length(obj.rgrid);
            obj.zindex=1:length(obj.zgrid);
            try
            obj.Rindex=h5partsquantity(obj.fullpath,hdf5group,'Rindex');
            obj.Zindex=h5partsquantity(obj.fullpath,hdf5group,'Zindex');
            
            obj.partindex=h5partsquantity(obj.fullpath,hdf5group,'partindex');
            catch
            end
            end
%             try
%                 obj.partindex = h5read(obj.fullpath,sprintf('%s/partindex',hdf5group));
%                 partindex=obj.partindex;
%                 partindex(partindex==-1)=NaN;
%                 %[~,Indices]=sort(partindex,'ascend');
%                 Indices=obj.partindex;
%                 for i=1:size(Indices,2)
%                     obj.H(Indices(:,i),i)=obj.H(:,i);
%                     obj.P(Indices(:,i),i)=obj.P(:,i);
%                     
%                     obj.R(Indices(:,i),i)=obj.R(:,i);
%                     obj.Z(Indices(:,i),i)=obj.Z(:,i);
%                     obj.Rindex(Indices(:,i),i)=obj.Rindex(:,i);
%                     obj.Zindex(Indices(:,i),i)=obj.Zindex(:,i);
%                     obj.VR(Indices(:,i),i)=obj.VR(:,i);
%                     obj.VZ(Indices(:,i),i)=obj.VZ(:,i);
%                     obj.VTHET(Indices(:,i),i)=obj.VTHET(:,i);
%                     obj.THET(Indices(:,i),i)=obj.THET(:,i);
%                 end
%             catch
%             end
%             clear partindex;
            
        end
        
        function quantity=H(obj,varargin)
            if(~iscell(varargin))
                indices=mat2cell(varargin);
            else
                indices=varargin;
            end
            if strcmp(indices{1},':')
                p=1:obj.VR.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            if size(indices,2)>2
                track=indices{3};
            else
                track=false;
            end
            quantity=0.5*obj.m*(obj.VR(p,t,track).^2+obj.VTHET(p,t,track).^2+obj.VZ(p,t,track).^2)+obj.partepot(p,t,track);
        end
        
        function quantity=P(obj,varargin)
            if(~iscell(varargin))
                indices=mat2cell(varargin);
            else
                indices=varargin;
            end
            if strcmp(indices{1},':')
                p=1:obj.R.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            if size(indices,2)>2
                track=indices{3};
            else
                track=false;
            end
            quantity=obj.R(p,t,track).*(obj.VTHET(p,t,track)*obj.m+obj.q*obj.parent.Atheta(obj.R(p,t,track),obj.Z(p,t,track)));
        end
        
        function quantity=Vpar(obj,varargin)
            %Vpar Computes the parallel velocity for the particle indices{1} at time indices{2}
            if(~iscell(varargin))
                indices=mat2cell(varargin);
            else
                indices=varargin;
            end
            if strcmp(indices{1},':')
                p=1:obj.R.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            if size(indices,2)>2
                track=indices{3};
            else
                track=false;
            end
            Zp=obj.Z(p,t,track);
            Rp=obj.R(p,t,track);
            Bzp=interp2(obj.zgrid,obj.rgrid,obj.parent.Bz',Zp,Rp);
            Brp=interp2(obj.zgrid,obj.rgrid,obj.parent.Br',Zp,Rp);
            Bp=interp2(obj.zgrid,obj.rgrid,obj.parent.B',Zp,Rp);
            costhet=Bzp./Bp;
            sinthet=Brp./Bp;
            quantity=obj.VR(p,t,track).*sinthet+obj.VZ(p,t,track).*costhet;
        end
        
        function quantity=Vperp(obj,varargin)
            %Vperp Computes the perpendicular velocity in the guidind center reference frame,
            % for the main specie particle indices{1} at time indices{2}
            
            if(~iscell(varargin))
                indices=mat2cell(varargin);
            else
                indices=varargin;
            end
            
            if strcmp(indices{1},':')
                p=1:obj.R.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            if size(indices,2)>2
                track=indices{3};
            else
                track=false;
            end
            if size(indices,2)>3
                gcs=indices{4};
            else
                gcs=false;
            end
            Zp=obj.Z(p,t,track);
            Rp=obj.R(p,t,track);
            Bzp=interp2(obj.zgrid,obj.rgrid,obj.parent.Bz',Zp,Rp);
            Brp=interp2(obj.zgrid,obj.rgrid,obj.parent.Br',Zp,Rp);
            Bp=interp2(obj.zgrid,obj.rgrid,obj.parent.B',Zp,Rp);
            costhet=Bzp./Bp;
            sinthet=Brp./Bp;
            Vdrift=zeros(size(Zp));
            if gcs
            for j=1:length(t)
                [~, tfield]=min(abs(obj.parent.t2d-obj.tpart(t(j))));
                timeEr=obj.parent.Er(:,:,tfield);
                timeEz=obj.parent.Ez(:,:,tfield);
                %posindE=sub2ind(size(timeEr),Rind(:,j),Zind(:,j));
                timeErp=interp2(obj.zgrid,obj.rgrid,timeEr,Zp(:,j),Rp(:,j));
                timeEzp=interp2(obj.zgrid,obj.rgrid,timeEz,Zp(:,j),Rp(:,j));
                Vdrift(:,j)=(timeEzp.*Brp(:,j)-timeErp.*Bzp(:,j))./Bp(:,j).^2;
            end
            end
            quantity=sqrt((obj.VTHET(p,t,track)-Vdrift).^2+(obj.VR(p,t,track).*costhet-obj.VZ(p,t,track).*sinthet).^2);
        end
    
        function sref = subsref(obj,s)
            % obj(i) is equivalent to obj.Data(i)
            switch s(1).type
                case '.'
                    if(strcmp(s(1).subs,'H'))
                            sref=H(obj,s(2).subs);
                    elseif(strcmp(s(1).subs,'P'))
                            sref=P(obj,s(2).subs);
                    else   
                        sref=builtin('subsref',obj,s);
                    end
                case '()'
                    sref=builtin('subsref',obj,s);
                case '{}'
                    error('MYDataClass:subsref',...
                        'Not a supported subscripted reference')
            end
        end
        
    end
end