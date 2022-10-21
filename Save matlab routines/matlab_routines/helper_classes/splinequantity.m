classdef splinequantity < h5quantity
    
    properties
        knotsr
        knotsz
        femorder
        index
        postscale
        rgrid
        zgrid
        invdr
        invdz
        rcenters
        zcenters
    end
    
    methods
        function obj=splinequantity(filename, dataset, knotsr, knotsz, femorder, scale, postscale, index)
            if nargin<6
                scale=1;
            end
            if nargin<8
                index=-1;
            end
            obj=obj@h5quantity(filename, dataset, length(knotsr)-2*femorder(2), length(knotsz)-2*femorder(1), scale);
            obj.knotsr=knotsr;
            obj.knotsz=knotsz;
            obj.femorder=double(femorder);
            if nargin < 7
                obj.postscale=ones(length(knotsr)-2*femorder(2), length(knotsz)-2*femorder(1));
            else
                obj.postscale=postscale;
            end
            if nargin < 8
                obj.index=-1;
            else
                obj.index=index;
            end
            obj.rgrid=obj.knotsr((1:obj.nr) +obj.femorder(2));
            obj.zgrid=obj.knotsz((1:obj.nz) +obj.femorder(1));
            [dz,dr]=meshgrid(obj.zgrid(3:end)-obj.zgrid(1:end-2),obj.rgrid(3:end)-obj.rgrid(1:end-2));
            obj.invdr=1./dr;
            obj.invdz=1./dz;
            obj.rcenters=movmean(obj.knotsr,femorder(2)+2);
            obj.rcenters=obj.rcenters(femorder(2):end-femorder(2));
            obj.zcenters=movmean(obj.knotsz,femorder(1)+2);
            obj.zcenters=obj.zcenters(femorder(1):end-femorder(1));
        end
        
        function quantity=coeffs(obj,indices)
            if strcmp(indices{1},':')
                r=1:obj.nr+obj.femorder(2)-1;
            else
                r=indices{1};
            end
            if strcmp(indices{2},':')
                z=1:obj.nz+obj.femorder(1)-1;
            else
                z=indices{2};
            end
            if strcmp(indices{3},':')
                t=1:obj.nt;
            else
                t=indices{3};
            end
            temp=zeros((obj.nz+obj.femorder(1)-1)*(obj.nr+obj.femorder(2)-1),length(t));
            if obj.index ~= -1
                if(length(unique(diff(t))) == 1 && length(t)>1)
                    stride=t(2)-t(1);
                    temp = h5read(obj.filename, obj.dataset,[obj.index 1 t(1)],[1 Inf length(t)],[1 1 stride]);
                else
                    for i=1:length(t)
                        temp(:,i) = h5read(obj.filename, obj.dataset,[obj.index 1 t(i)],[1 Inf 1]);
                    end
                end
            else
                if(sum(diff(diff(t))) == 0 && length(t)>1)
                    stride=t(2)-t(1);
                    temp = h5read(obj.filename, obj.dataset,[1 t(1)],[Inf length(t)],[1 stride]);
                else
                    for i=1:length(t)
                        temp(:,i) = h5read(obj.filename, obj.dataset,[1 t(i)],[Inf 1]) ;
                    end
                end
            end
            temp=reshape(squeeze(temp),obj.nz+obj.femorder(1)-1,obj.nr+obj.femorder(2)-1,[]);
            temp=permute(temp,[2,1,3]);
            quantity=temp(r,z,:)*obj.scale;
        end
        
        function quantity=val(obj,indices)
            if strcmp(indices{1},':')
                r=1:obj.nr;
            else
                r=indices{1};
            end
            if strcmp(indices{2},':')
                z=1:obj.nz;
            else
                z=indices{2};
            end
            if strcmp(indices{3},':')
                t=1:obj.nt;
            else
                t=indices{3};
            end
            count=length(t);
            temp=obj.coeffs({':',':',t});
            quantity=zeros(length(r),length(z),count);
            [Z,R]=meshgrid(obj.zcenters,obj.rcenters);
            [zg,rg]=meshgrid(obj.zgrid,obj.rgrid);
                for i=1:size(temp,3)
                    quantity(:,:,i)=obj.postscale(r,z).*fnval(spmak({obj.knotsr,obj.knotsz},temp(:,:,i)),{obj.knotsr(r+obj.femorder(2)),obj.knotsz(z+obj.femorder(1))});
                    %valued=interp2(Z,R,squeeze(temp(:,:,i)),zg,rg);
                    %quantity(:,:,i)=valued(r,z);
                end      
        end
        
        function quantity=posval(obj,indices)
            if strcmp(indices{1},':')
                r=1:obj.nr;
            else
                r=indices{1};
            end
            if strcmp(indices{2},':')
                z=1:obj.nz;
            else
                z=indices{2};
            end
            if strcmp(indices{3},':')
                t=1:obj.nt;
            else
                t=indices{3};
            end
            count=length(t);
            temp=obj.coeffs({':',':',t});
            if(length(r) ~= length(z))
                error("r and z array must be the same size")
            end
            quantity=zeros(min(length(r),length(z)),count);
            for i=1:size(temp,3)
                quantity(:,i)=fnval(spmak({obj.knotsr,obj.knotsz},temp(:,:,i)),[r(:)';z(:)']);
            end
                
        end
        
        
        function quantity=der(obj,indices)
            if strcmp(indices{1},':')
                r=1:obj.nr;
            else
                r=indices{1};
            end
            if strcmp(indices{2},':')
                z=1:obj.nz;
            else
                z=indices{2};
            end
            if strcmp(indices{3},':')
                t=1:obj.nt;
            else
                t=indices{3};
            end
            order=indices{4};
            count=length(t);
            temp=obj.coeffs({':',':',t});
            quantity=zeros(length(r),length(z),count);
            for i=1:size(temp,3)
                %f=spmak({obj.knotsr,obj.knotsz},temp(:,:,i));
                %preder=fnval(f,{obj.knotsr(r+obj.femorder(2)),obj.knotsz(z+obj.femorder(1))});
                preder=obj.val({r,z,t(i)});%
                    if order(1)>0
                        preder(2:end-1,2:end-1)=(preder(3:end,2:end-1)-preder(1:end-2,2:end-1)).*obj.invdr;
                    end
                    if order(2)>0
                        preder(2:end-1,2:end-1)=(preder(2:end-1,3:end)-preder(2:end-1,1:end-2)).*obj.invdz;
                    end
                    quantity(2:end-1,2:end-1,i)=preder(2:end-1,2:end-1);
            end
        end
    end
end
