classdef h5partsquantity
    
    properties
        filename
        group
        dataset
        nparts
        nt
        scale
    end
    
    methods
        function obj=h5partsquantity(filename, group, dataset, scale)
            obj.filename=filename;
            obj.dataset=dataset;
            obj.group=group;
            obj.nparts=h5info(filename,sprintf('%s/%s',group,dataset)).Dataspace.Size(1);
            obj.nt=h5info(filename,sprintf('%s/%s',group,dataset)).Dataspace.Size(end);
            if nargin < 4
                obj.scale=1;
            else
                obj.scale=scale;
            end
        end
        
        function ind=end(obj,k,n)
            switch k
                case 1
                    ind=obj.nparts;
                case 2
                    ind=obj.nt;
                case default
                    error('Invalid number of dimensions');
            end
        end
        
        function sref = subsref(obj,s)
            % obj(i) is equivalent to obj.Data(i)
            switch s(1).type
                case '.'
                    if(strcmp(s(1).subs,'val'))
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
        
         function quantity=val(obj,indices)
            if strcmp(indices{1},':')
                p=1:obj.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:obj.nt;
            else
                t=indices{2};
            end
            if size(indices,2)>2
                track=indices{3};
            else
                track=false;
            end
            
            if track
                quantity=sparse(max(p),length(t));
            for i=1:length(t)
                 indices = h5read(obj.filename, sprintf('%s/%s',obj.group,'partindex'),[1 t(i)],[Inf 1]);
                 indices(indices<=0)=max(indices)+1;
                 quantity(indices,i) = h5read(obj.filename, sprintf('%s/%s',obj.group,obj.dataset),[1 t(i)],[Inf 1]); 
            end
            else
                quantity=zeros(length(p),length(t));
                for i=1:length(t)
                 quantity(:,i) = h5read(obj.filename, sprintf('%s/%s',obj.group,obj.dataset),[1 t(i)],[length(p) 1]); 
                end
            end
            quantity=quantity(p,:)*obj.scale;
        end
    end
end
