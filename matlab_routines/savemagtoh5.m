function savemagtoh5(filename,r,z,Athet,Br,Bz,overwrite)
% saves a magnetic field to be used with Fennecs
% if Br or Bz is not set it is calulated from the magnetic vector potential
% Overwrite defines if the .h5 file must be overwritten
    if nargin<6
        [Bz,~]=-gradient(Athet,mean(diff(z)));
    end
    if nargin<5
        [~, Br]=gradient(r.*Athet,1,mean(diff(r)));
        Br=Br./r;
    end 
    if nargin<7
        overwrite=false;
    end
    
    if overwrite && exist(filename,'file')==2
        delete(filename);
    end
    
    h5create(filename,'/mag/r',length(r));
    h5create(filename,'/mag/z',length(z));
    h5create(filename,'/mag/Athet',size(Athet'));
    h5create(filename,'/mag/Br',size(Br'));
    h5create(filename,'/mag/Bz',size(Bz'));
    h5write(filename,'/mag/r',r) % stores the radial grid position
    h5write(filename,'/mag/z',z) % stores the axial grid position
    h5write(filename,'/mag/Athet',Athet')
    h5write(filename,'/mag/Br',Br')
    h5write(filename,'/mag/Bz',Bz')
    
end