classdef h5quantity
    
    properties
        filename
        dataset
        nr
        nz
        nt
        scale
    end
    
    methods
        function obj=h5quantity(filename, dataset, nr, nz, scale)
            obj.filename=filename;
            obj.nr=nr;
            obj.nz=nz;
            obj.dataset=dataset;
            obj.nt=h5info(filename,dataset).Dataspace.Size(end);
            if nargin < 5
                obj.scale=1;
            else
                obj.scale=scale;
            end
        end
        
        function ind=end(obj,k,n)
            switch k
                case 1
                    ind=obj.nr;
                case 2
                    ind=obj.nz;
                case 3
                    ind=obj.nt;
                case default
                    error('Invalid number of dimensions');
            end
        end
        
        function sref = subsref(obj,s)
            % obj(i) is equivalent to obj.Data(i)
            switch s(1).type
                case '.'
                    if(strcmp(s(1).subs,'coeffs') || strcmp(s(1).subs,'val')|| strcmp(s(1).subs,'der'))
                        sref = obj.(s(1).subs)(s(2).subs);
                    else
                        sref=builtin('subsref',obj,s);
                    end
                case '()'
                    sref = val(obj,s(1).subs);
                case '{}'
                    error('MYDataClass:subsref',...
                        'Not a supported subscripted reference')
            end
        end
    end
end
