classdef splinedensity < splinequantity
    
    properties
        invVolume
    end
    
    methods
        function obj=splinedensity(filename, dataset, knotsr, knotsz, femorder, Volume, MacropartSize, postscale, index)
            obj=obj@splinequantity(filename, dataset, knotsr, knotsz, femorder, MacropartSize, postscale, index);
            obj.invVolume=1./Volume;
            obj.invVolume(isinf(obj.invVolume))=0;
        end
        
        function quantity=coeffs(obj,indices)
            quantity=coeffs@splinequantity(obj,indices);
            for i=1:size(quantity,3)
                quantity(:,:,i)=quantity(:,:,i).*obj.invVolume;
            end
        end
    end
end