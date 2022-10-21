classdef splineenergy < h5quantity
    
    properties
        invVolume
        knotsr
        knotsz
        femorder
        postscale
    end
    
    methods (Access=public)
        function obj=splineenergy(filename, dataset, knotsr, knotsz, femorder, Volume, vnorm, postscale)
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
        end
        
        function quantity=coeffs(obj,indices)
            if strcmp(indices{1},':')
                ij=1:3;
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
                    quantity(tempi,:,:,timei)=permute(reshape(temp(:,timei),obj.nz+obj.femorder(1)-1,obj.nr+obj.femorder(2)-1),[2,1,3]);
                end
            end
            quantity=quantity*obj.scale;
        end
        
        function quantity=val(obj,indices)
            if strcmp(indices{1},':')
                ij=1:3;
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
            quantity=zeros(length(ij),length(r),length(z),count);
            for j=1:size(temp,1)
                for i=1:size(temp,4)
                    vali=obj.postscale.*fnval(spmak({obj.knotsr,obj.knotsz},squeeze(temp(j,:,:,i))),{obj.knotsr(r+obj.femorder(2)),obj.knotsz(z+obj.femorder(1))});
                    quantity(j,:,:,i)=vali(r,z);
                end
            end
        end
        
        function quantity=der(obj,indices)
            if strcmp(indices{1},':')
                ij=1:3;
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
                    quantity(j,:,:,i)=fnval(fnder(...
                        spmak({obj.knotsr,obj.knotsz},squeeze(temp(j,:,:,i)))...
                        ,order),{obj.knotsr(r+obj.femorder(2)),obj.knotsz(z+obj.femorder(1))});
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
        function quantity = readtensorindex(obj, ii, t)
            if strcmp(ii,':') || length(ii)>1
                error('Unable to read several indices at the same time');
            end
            if strcmp(t,':')
                t=1:obj.nt;
            end
            vii=zeros((obj.nz+obj.femorder(1)-1)*(obj.nr+obj.femorder(2)-1),length(t));
            n=vii;
            
            switch ii
                case 1
                    id=5;
                case 2
                    id=8;
                case 3
                    id=10;
                case default
                    error('Invalid Energy index'); 
            end
           
            for k=1:length(t)
                vii(:,k) = h5read(obj.filename, obj.dataset,[id   1 t(k)],[1 Inf 1]);
                n(:,k)   = h5read(obj.filename, obj.dataset,[1    1 t(k)],[1 Inf 1]);
            end
            n=1./n;
            n(isinf(n))=0;
            quantity=vii.*n;
        end
    end
end