classdef splinepressure < h5quantity
    
    properties
        invVolume
        knotsr
        knotsz
        femorder
        postscale
        rgrid
        zgrid
        invdr
        invdz
        rcenters
        zcenters
    end
    
    methods (Access=public)
        function obj=splinepressure(filename, dataset, knotsr, knotsz, femorder, Volume, vnorm, postscale)
            if nargin<7
                vnorm=1;
            end
            obj=obj@h5quantity(filename, dataset, length(knotsr)-2*femorder(2), length(knotsz)-2*femorder(1), vnorm);
            obj.knotsr=knotsr;
            obj.knotsz=knotsz;
            obj.femorder=double(femorder);
            obj.invVolume=1./Volume;
            obj.invVolume(isinf(obj.invVolume))=0;
            if nargin <8
                obj.postscale=1;
            else
                obj.postscale=postscale;
            end
            obj.rgrid=obj.knotsr((1:obj.nr) +obj.femorder(2));
            obj.zgrid=obj.knotsz((1:obj.nz) +obj.femorder(1));
            [dr,dz]=meshgrid(obj.rgrid(3:end)-obj.rgrid(1:end-2),obj.zgrid(3:end)-obj.zgrid(1:end-2));
            obj.invdr=1./dr';
            obj.invdz=1./dz';
            obj.rcenters=movmean(obj.knotsr,femorder(2)+2);
            obj.rcenters=obj.rcenters(femorder(2):end-femorder(2));
            obj.zcenters=movmean(obj.knotsz,femorder(1)+2);
            obj.zcenters=obj.zcenters(femorder(1):end-femorder(1));
        end
        
        function quantity=coeffs(obj,indices)
            if strcmp(indices{1},':')
                ij=1:6;
            else
                ij=indices{1};
            end
            if strcmp(indices{2},':')
                r=1:obj.nr+obj.femorder(2)-1;
            else
                r=indices{2};
            end
            if strcmp(indices{3},':')
                z=1:obj.nz+obj.femorder(1)-1;
            else
                z=indices{3};
            end
            if strcmp(indices{4},':')
                t=1:obj.nt;
            else
                t=indices{4};
            end
            quantity=zeros(length(ij),(obj.nr+obj.femorder(2)-1),(obj.nz+obj.femorder(1)-1),length(t));
            %temp=zeros((obj.nz+obj.femorder(1)-1)*(obj.nr+obj.femorder(2)-1),length(t));
            for tempi=1:length(ij)
                temp=obj.readtensorindex(ij(tempi),t);
                for timei=1:size(temp,2)
                    quantity(tempi,:,:,timei)=permute(reshape(temp(:,timei),obj.nz+obj.femorder(1)-1,obj.nr+obj.femorder(2)-1),[2,1,3]).*obj.invVolume;
                end
            end
            quantity=quantity*obj.scale;
        end
        
        function quantity=val(obj,indices)
            if strcmp(indices{1},':')
                ij=1:6;
            else
                ij=indices{1};
            end
            if strcmp(indices{2},':')
                r=1:obj.nr;
            else
                r=indices{2};
            end
            if strcmp(indices{3},':')
                z=1:obj.nz;
            else
                z=indices{3};
            end
            if strcmp(indices{4},':')
                t=1:obj.nt;
            else
                t=indices{4};
            end
            count=length(t);
            temp=obj.coeffs({ij,':',':',t});
            quantity=zeros(length(ij),max(r)-min(r)+1,max(z)-min(z)+1,count);
            %[Z,R]=meshgrid(obj.zcenters,obj.rcenters);
            %[zg,rg]=meshgrid(obj.zgrid,obj.rgrid);
            for j=1:size(temp,1)
                for i=1:size(temp,4)
                    quantity(j,:,:,i)=fnval(spmak({obj.knotsr,obj.knotsz},squeeze(temp(j,:,:,i))),{obj.knotsr(r+obj.femorder(2)),obj.knotsz(z+obj.femorder(1))});
                    %quantity(j,:,:,i)=interp2(Z,R,squeeze(temp(j,:,:,i)),zg,rg);
                end
            end
        end
        
        function quantity=der(obj,indices)
            if strcmp(indices{1},':')
                ij=1:6;
            else
                ij=indices{1};
            end
            if strcmp(indices{2},':')
                r=1:obj.nr;
            else
                r=indices{2};
            end
            if strcmp(indices{3},':')
                z=1:obj.nz;
            else
                z=indices{3};
            end
            if strcmp(indices{4},':')
                t=1:obj.nt;
            else
                t=indices{4};
            end
            if length(indices)<5
                order=[1,1];
            else
                order=indices{5};
            end
            count=length(t);
            temp=obj.coeffs({ij,':',':',t});
            quantity=zeros(length(ij),max(r)-min(r)+1,max(z)-min(z)+1,count);
            for j=1:size(temp,1)
                for i=1:size(temp,4)
                    %preder=fnval(...
                    %           spmak({obj.knotsr,obj.knotsz},squeeze(temp(j,:,:,i)))...
                    %    ,{obj.knotsr(r+obj.femorder(2)),obj.knotsz(z+obj.femorder(1))});
                    preder=squeeze(obj.val({ij(j),r,z,t(i)}));
                    if order(1)>0
                        preder(2:end-1,2:end-1)=(preder(3:end,2:end-1)-preder(1:end-2,2:end-1)).*obj.invdr;
                    end
                    if order(2)>0
                        preder(2:end-1,2:end-1)=(preder(2:end-1,3:end)-preder(2:end-1,1:end-2)).*obj.invdz;
                    end
                    quantity(j,2:end-1,2:end-1,i)=preder(2:end-1,2:end-1);
                end
            end
        end
        
        function ind=end(obj,k,n)
            switch k
                case 1
                    ind=6;
                case 2
                    ind=obj.nr;
                case 3
                    ind=obj.nz;
                case 4
                    ind=obj.nt;
                case default
                    error('Invalid number of dimensions');
            end
        end
    end
    
    methods (Access=private)
        function quantity = readtensorindex(obj, ij, t)
            if strcmp(ij,':') || length(ij)>1
                error('Unable to read several indices at the same time');
            end
            if strcmp(t,':')
                t=1:obj.nt;
            end
            ui=zeros((obj.nz+obj.femorder(1)-1)*(obj.nr+obj.femorder(2)-1),length(t));
            uj=ui;
            n=ui;
            vij=ui;
            
            switch ij
                case 1
                    i=1;
                    j=1;
                case 2
                    i=1;
                    j=2;
                case 3
                    i=1;
                    j=3;
                case 4
                    i=2;
                    j=2;
                case 5
                    i=2;
                    j=3;
                case 6
                    i=3;
                    j=3;
                case default
                    error('Invalid Pressure index'); 
            end
           
            for k=1:length(t)
                vij(:,k) = h5read(obj.filename, obj.dataset,[ij+4 1 t(k)],[1 Inf 1]);
                ui(:,k)  = h5read(obj.filename, obj.dataset,[i+1  1 t(k)],[1 Inf 1]);
                uj(:,k)  = h5read(obj.filename, obj.dataset,[j+1  1 t(k)],[1 Inf 1]);
                n(:,k)   = h5read(obj.filename, obj.dataset,[1    1 t(k)],[1 Inf 1]);
            end
            n=1./n;
            n(isinf(n))=0;
            nv2=vij;
            nuiuj=ui.*uj.*n;
            quantity=squeeze(nv2-nuiuj);
        end
    end
end