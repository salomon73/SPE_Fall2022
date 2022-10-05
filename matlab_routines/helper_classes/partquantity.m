classdef partquantity < h5quantity
    
    properties
        
    end
    
    methods
        function obj=partquantity(filename, dataset,scale)
            if nargin<3
                scale=1;
            end
            obj=obj@h5quantity(filename, dataset, 0, 0, scale);
            
        end
        
        function quantity=val(obj,indices)
            if indices{1}==':'
                r=1:obj.nr;
            else
                r=indices{1};
            end
            if indices{2}==':'
                z=1:obj.nz;
            else
                z=indices{2};
            end
            if indices{3}==':'
                t=1:obj.nt;
            else
                t=indices{3};
            end
            if obj.centered
                 quantity=zeros((obj.nz-1)*(obj.nr-1),length(t));
            else
                quantity=zeros((obj.nz)*(obj.nr),length(t));
            end
            for i=1:length(t)
                 quantity(:,i) = h5read(obj.filename, obj.dataset,[1 t(i)],[Inf 1]) ;
            end
            if obj.centered
                 quantity=reshape(quantity,obj.nz-1,obj.nr-1,[]);
            else
                quantity=reshape(quantity,obj.nz,obj.nr,[]);
            end
            quantity=permute(quantity,[2,1,3]);
            
            if obj.centered
                quantity=padarray(quantity,[1,1],0,'post');
            end
            quantity=quantity(r,z,:)*obj.scale;
        end
    end
end
