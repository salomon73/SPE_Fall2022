function savegeomtoh5(filename,splines,dist_extent,overwrite)
% saves the geometry to be used with Fennecs
% splines is a structure array containing each individual boundary
% Overwrite defines if the .h5 file must be overwritten
% dist_extent is the distance in m on which the geometric weight goes from 0 to 1 
% dist_extent can also be set in the job input file
% To look how to define each spline read the function
if nargin<4
        overwrite=false;
end    
if nargin<4
        epsge=1e-6;
    end
    group='/geometry_spl/';
    try
        fid=H5F.create(filename);
        H5F.flush(fid,'H5F_SCOPE_LOCAL')
    catch
        if overwrite
            delete(filename);
            fid=H5F.create(filename);
            H5F.flush(fid,'H5F_SCOPE_LOCAL')
        else
            error('cFile "%s" already exists\n',filename) ;
        end
    end
    nbsplines=length(splines);
    gid=H5G.create(fid,group,64);
    h5writeatt(filename,group,'nbsplines',nbsplines); % total number of spline curves
    h5writeatt(filename,group,'dist_extent',dist_extent);
    if iscell(splines)
        splines=cell2mat(splines);
    end
    for i=1:nbsplines
        grp=sprintf('/geometry_spl/%2.2d',i);
        gid=H5G.create(fid,grp,64);
        h5writeatt(filename,grp,'Dirichlet_val',splines(i).Dval); % value of the electic potential in V at the boundary surface
        h5writeatt(filename,grp,'order',splines(i).order); % order of the spline used to define the boundary
        h5writeatt(filename,grp,'dim',splines(i).dim); % dimension of the spline curve
        % epsge and epsce have no clear effect default values can be used
        h5writeatt(filename,grp,'epsge',splines(i).epsge); % geometric presicion used by SISL for calculating distance to curve
        h5writeatt(filename,grp,'epsce',splines(i).epsce); % numeric presicion used by SISL for calculating distance to curve
        h5writeatt(filename,grp,'name',splines(i).name); % name of the boundary for readability of h5(not used in fennecs)
        h5writeatt(filename,grp,'type',splines(i).type); % boundary condition type (0 dirichlet fixed value)
                                                         % 1 Dirichlet with value depending  on s
                                                         % 2 Natural
        h5writeatt(filename,grp,'periodic',splines(i).periodic); % periodicity of the curve 1 open, 0 closed, -1 closed periodic
        points=splines(i).points;
        if size(splines(i).points,2)>splines(i).dim
            points=points';
        end
        h5create(filename,[grp,'/pos'],size(points));
        h5write(filename,[grp,'/pos'],points); % control points defining the curve
    end
    
    H5F.close(fid);
end