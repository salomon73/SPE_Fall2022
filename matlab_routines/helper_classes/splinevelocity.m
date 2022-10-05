classdef splinevelocity < splinequantity
    
    properties
        Partspercell
    end
    
    methods
        function obj=splinevelocity(filename, dataset, knotsr, knotsz, femorder, vnorm, postscale, index)
            obj@splinequantity(filename, dataset, knotsr, knotsz, femorder, vnorm, postscale, index);
            obj.Partspercell=splinequantity(filename, dataset, knotsr, knotsz, femorder, 1, 1, 1);
        end
        
        function quantity=coeffs(obj,indices)
            scalefactgrid=coeffs@splinequantity(obj.Partspercell,indices);
            scalefactgrid=1./scalefactgrid;
            scalefactgrid(isinf(scalefactgrid))=0;
            % We divide by the local N top obtain the velocity
            quantity=coeffs@splinequantity(obj,indices).*scalefactgrid;
        end
    end
end