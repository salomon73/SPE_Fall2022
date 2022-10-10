classdef griddensity < gridquantity
    
    properties
        invVolume
    end
    
    methods
        function obj=griddensity(filename, dataset, nr, nz, Volume, MacropartSize, centered)
            if nargin<6
                MacropartSize=1;
            end
            if nargin<7
                centered=false;
            end
            obj=obj@gridquantity(filename, dataset, nr, nz, MacropartSize, centered);
            obj.invVolume=1./Volume;
            obj.invVolume(isinf(obj.invVolume))=0;
            if centered
                obj.invVolume=padarray(obj.invVolume,[1,1],0,'post');
            end
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
            quantity=val@gridquantity(obj,indices);
            for i=1:size(quantity,3)
                quantity(:,:,i)=quantity(:,:,i).*obj.invVolume(r,z);
            end
        end
    end
end