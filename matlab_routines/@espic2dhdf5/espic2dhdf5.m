classdef espic2dhdf5
    %espic2dhdf5 General class used to treat hdf5 result files of espic2d code
    %   A result file is loaded with a call to M=espic2dhdf5(filename) where filename is the relative or absolute file path
    %   after loading, several quantities and composite diagnostics such as moments of the distribution function or individual particles
    %   quantities can be accessed.
    properties
        filename
        name
        folder
        fullpath
        timestamp
        info
        t0d
        t1d
        t2d
        tpart
        it0
        it1
        it2
        restartsteps
        restarttimes
        
        %% Physical constants
        vlight=299792458;
        qe=1.60217662E-19;
        me=9.109383E-31;
        eps_0=8.85418781762E-12;
        kb=1.38064852E-23;
        
        %% Run parameters
        dt    % simulation time step
        nrun  % number of time steps simulated
        nlres
        nlsave
        nlclassical % Was the equation of motion solved in the classical framework
        nlPhis % Was the self-consistent electric field computed
        nz % number of intervals in the z direction for the grid
        nnr % number of intervals in the r direction for the grid for each of the 3 mesh regions
        lz % physical axial dimension of the simulation space
        nplasma % Number of initial macro particles
        potinn % Normalized electric potential at the coaxial insert
        potout % Normalized electric potential at the cylinder surface
        B0     % Normalization for the magnetic field
        Rcurv  % Magnetic mirror ratio
        width  % Magnetic mirror length
        n0     % Initial particle density in case of old particle loading
        temp   % Initial particle temperature in case of old particle loading
        femorder % finite element method order in z and r direction
        ngauss   % Order of the Gauss integration method for the FEM
        plasmadim % initial dimensions of the plasma for the old particle loading system
        radii     % Radial limits of the three mesh regions coarse,fine,coarse
        H0        % Initial particle Energy for Davidsons distribution function
        P0        % Initial particle Angular momentum for Davidsons distribution function
        normalized % Are the parts quantities normalized in the h5 file
        nbspecies  % Number of species simulated
        
        
        %% Frequencies
        omepe  % Reference plasma frequency used for normalization
        omece  % Reference cyclotronic frequency for normalization
        
        %% Normalizations
        tnorm    % Time normalization
        rnorm    % Dimension normalization
        bnorm    % Magnetic field normalization
        enorm    % Electric field normalization
        phinorm  % Electric potential normalization
        vnorm    % Velocity normalization
        
        %% Grid data
        rgrid % Radial grid position points
        zgrid % Axial grid position points
        dz    % Axial grid step
        dr    % Radial grid step for the three mesh regions
        CellVol % Volume of the cell used for density calculation
        celltype % type of cell -1 outside 1 inside 0 border
        linked_s % location of linked spline
        bsplinetype
        
        %% Magnetic field
        Br % Radial magnetic field
        Bz % Axial magnetic field
        Athet % Azimuthal component of the Magnetic potential vector
        rAthet % r*Athet used for the representation of magnetic field lines
        B      % Magnetic field amplitude
        sinthet % ratio to project quantities along the magnetic field lines
        costhet % ratio to project quantities along the magnetic field lines
        
        %% Energies
        epot  % Time evolution of the particles potential energy
        ekin  % Time evolution of the particles kinetic energy
        etot  % Time evolution of the particles total energy
        etot0 % Time evolution of the reference particle total energy
        eerr  % Time evolution of the error on the energy conservation
        npart % Time evolution of the number of simulated
        
        %% 2D time data evaluated on grid points
        N         % main specie Density
        fluidUR   % main specie radial fluid velocity
        fluidUZ   % main specie axial fluid velocity
        fluidUTHET % main specie azimuthal fluid velocity
        pot        % Electric potential evaluated at grid points
        potxt     % External Electric potential evaluated at grid points
        phi        % Electric potential in spline form
        Er         % Radial electric field
        Ez         % Axial electric field
        Erxt      % External Radial electric field
        Ezxt      % External Axial electric field
        Presstens  % Pressure tensor
        fluidEkin  % average kinetic energy in each direction
        
        %% Splines
        knotsr     % Spline radial knots
        knotsz     % Spline axial knots
        
        %% Particle parameters
        weight  % Macro particle numerical weight of the main specie
        qsim    % Macro particle charge
        msim    % Macro particle mass
        nbparts   % Time evolution of the number of simulated particles
        partepot % Electric potential at the particles positions
        R        % Particles radial position
        Z        % Particles axial position
        Rindex   % Particles radial grid index
        Zindex   % Particles axial grid index
        partindex % Particles unique id for tracing trajectories
        VR        % Particles radial velocity
        VZ        % Particles axial velocity
        VTHET     % Particles azimuthal velocity
        THET      % Particles azimuthal position
        species   % Array containing the other simulated species
        
        %% Celldiag
        celldiag  % Array containing the cell diagnostic data
        nbcelldiag % Total number of cell diagnostics
        
        %% Curvilinear geometry
        conformgeom % stores if we use the conforming or nonconforming boundary conditions
        r_a
        r_b
        z_r
        z_0
        r_0
        r_r
        L_r
        L_z
        Interior
        above1
        above2
        interior
        walltype
        geomweight
        dirichletweight
        gtilde
        spl_bound
        
        
        
        %% Maxwell source parameters
        maxwellsrce
        
        %% Collision with neutral parameters
        neutcol
        nudcol % effective momentum collision frequency
        
        %% Non ideal power supply
        psupply
        
    end
    
    methods
        function file=file(obj)
            % returns the h5 file name
            file=obj.filename;
        end
        
        function obj = espic2dhdf5(filename,readparts,old)
            % Reads the new result file filename and read the parts data if readparts==true
            
            % adds the helper_classes folder to the path
            matlabfuncpath = dir([mfilename('fullpath'),'.m']);
            addpath(sprintf('%s/../helper_classes',matlabfuncpath.folder));
            addpath(sprintf('%s/../extrema',matlabfuncpath.folder));
            addpath(sprintf('%s/../export_fig',matlabfuncpath.folder));
            rehash path
            % Try catch are there for compatibility with older simulation files
            filedata=dir(filename);
            if (isempty(filedata))
                error("File: ""%s"" doesn't exist",filename)
            end
            obj.folder=filedata.folder;
            obj.filename=filename;
            [~, obj.name, ext] = fileparts(obj.filename);
            obj.filename=[obj.name,ext];
            obj.fullpath=[obj.folder,'/',obj.filename];
            obj.timestamp=filedata.date;
            if nargin==1
                readparts=true;
            end
            if nargin<3
                old=false;
            end
            %obj.info=h5info(filename);
            
            
            %% Read the run parameters
            obj.dt = h5readatt(obj.fullpath,'/data/input.00/','dt');
            obj.nrun = h5readatt(obj.fullpath,'/data/input.00/','nrun');
            obj.nlres = strcmp(h5readatt(obj.fullpath,'/data/input.00/','nlres'),'y');
            obj.nlsave = strcmp(h5readatt(obj.fullpath,'/data/input.00/','nlsave'),'y');
            obj.nlclassical =strcmp(h5readatt(obj.fullpath,'/data/input.00/','nlclassical'),'y');
            obj.nlPhis =strcmp(h5readatt(obj.fullpath,'/data/input.00/','nlPhis'),'y');
            obj.nz = h5readatt(obj.fullpath,'/data/input.00/','nz');
            obj.nnr = h5read(obj.fullpath,'/data/input.00/nnr');
            obj.lz = h5read(obj.fullpath,'/data/input.00/lz');
            obj.qsim = h5readatt(obj.fullpath,'/data/input.00/','qsim');
            obj.msim = h5readatt(obj.fullpath,'/data/input.00/','msim');
            
            try
                obj.r_a=h5readatt(obj.fullpath,'/data/input.00/geometry','r_a');
                obj.r_b=h5readatt(obj.fullpath,'/data/input.00/geometry','r_b');
                obj.z_r=h5readatt(obj.fullpath,'/data/input.00/geometry','z_r');
                obj.r_r=h5readatt(obj.fullpath,'/data/input.00/geometry','r_r');
                obj.r_0=h5readatt(obj.fullpath,'/data/input.00/geometry','r_0');
                obj.z_0=h5readatt(obj.fullpath,'/data/input.00/geometry','z_0');
                obj.above1=h5readatt(obj.fullpath,'/data/input.00/geometry','above1');
                obj.above2=h5readatt(obj.fullpath,'/data/input.00/geometry','above2');
                obj.interior=h5readatt(obj.fullpath,'/data/input.00/geometry','interior');
                obj.walltype=h5readatt(obj.fullpath,'/data/input.00/geometry','walltype');
                try
                    obj.L_r=h5readatt(obj.fullpath,'/data/input.00/geometry','L_r');
                    obj.L_z=h5readatt(obj.fullpath,'/data/input.00/geometry','L_z');
                catch
                end
                obj.conformgeom=false;
            catch
                obj.conformgeom=true;
                obj.walltype=0;
                obj.r_a=obj.rgrid(1);
                obj.r_b=obj.rgrid(end);
                obj.above1=1;
                obj.above2=-1;
                obj.L_r=0;
                obj.L_z=0;
            end
            
            try
                obj.weight=h5readatt(obj.fullpath,'/data/part/','weight');
            catch
                obj.weight=obj.msim/obj.me;
            end
            filesgrpinfo=h5info(obj.fullpath,'/files');
            nbrst=h5readatt(obj.fullpath,'/files','jobnum');
            obj.restartsteps(1)=0;
            obj.restarttimes(1)=0;
            grp=sprintf('/data/input.%02i/',0);
            obj.dt(1)=h5readatt(obj.fullpath,grp,'dt');
            for i=1:nbrst
               grp=sprintf('/data/input.%02i/',i);
               obj.restartsteps(i+1)= h5readatt(obj.fullpath,grp,'startstep');
               obj.restarttimes(i+1)= obj.restarttimes(i) + (obj.restartsteps(i+1)-obj.restartsteps(i))*obj.dt(i);
               obj.dt(i+1)=h5readatt(obj.fullpath,grp,'dt');
            end
            obj.nplasma = h5readatt(obj.fullpath,'/data/input.00/','nplasma');
            obj.potinn = h5readatt(obj.fullpath,'/data/input.00/','potinn');
            obj.potout = h5readatt(obj.fullpath,'/data/input.00/','potout');
            obj.B0 = h5readatt(obj.fullpath,'/data/input.00/','B0');
            obj.Rcurv = h5readatt(obj.fullpath,'/data/input.00/','Rcurv');
            obj.width = h5readatt(obj.fullpath,'/data/input.00/','width');
            obj.n0 = h5readatt(obj.fullpath,'/data/input.00/','n0');
            obj.temp = h5readatt(obj.fullpath,'/data/input.00/','temp');
            try
                obj.it0 = h5readatt(obj.fullpath,'/data/input.00/','it0d');
                obj.it1 = h5readatt(obj.fullpath,'/data/input.00/','it2d');
                obj.it2 = h5readatt(obj.fullpath,'/data/input.00/','itparts');
            catch
                obj.it0 = h5readatt(obj.fullpath,'/data/input.00/','it0');
                obj.it1 = h5readatt(obj.fullpath,'/data/input.00/','it1');
                obj.it1 = h5readatt(obj.fullpath,'/data/input.00/','it2');
            end
            try
                try
                    obj.nbspecies=h5readatt(obj.fullpath,'/data/part/','nbspecies');
                catch
                    obj.nbspecies=h5readatt(obj.fullpath,'/data/input.00/','nbspecies');
                end
                obj.normalized=strcmp(h5readatt(obj.fullpath,'/data/input.00/','rawparts'),'y');
            catch
                obj.nbspecies=1;
                obj.normalized=false;
            end
            try
                obj.nbcelldiag=h5readatt(obj.fullpath,'/data/celldiag/','nbcelldiag');
            catch
                obj.nbcelldiag=0;
            end
            
            obj.omepe=sqrt(abs(obj.n0)*obj.qe^2/(obj.me*obj.eps_0));
            obj.omece=obj.qe*obj.B0/obj.me;
            
            obj.npart= h5read(obj.fullpath, '/data/var0d/nbparts');
            try
                obj.nudcol= h5read(obj.fullpath, '/data/var0d/nudcol');
            catch
            end
            try
                obj.H0 = h5read(obj.fullpath,'/data/input.00/H0');
                obj.P0 = h5read(obj.fullpath,'/data/input.00/P0');
            catch
                obj.H0=3.2e-14;
                obj.P0=8.66e-25;
            end
            
            % Normalizations
            if old
                obj.tnorm=abs(1/obj.omepe);
            else
                obj.tnorm=min(abs(1/obj.omepe),abs(1/obj.omece));
            end
            obj.rnorm=obj.vlight*obj.tnorm;
            obj.bnorm=obj.B0;
            obj.enorm=obj.vlight*obj.bnorm;
            obj.phinorm=obj.enorm*obj.rnorm;
            obj.vnorm=obj.vlight;
            
            % Grid data
            obj.rgrid= h5read(obj.fullpath, '/data/var1d/rgrid')*obj.rnorm;
            obj.zgrid= h5read(obj.fullpath, '/data/var1d/zgrid')*obj.rnorm;
            obj.dz=(obj.zgrid(end)-obj.zgrid(1))/double(obj.nz);
            rid=1;
            for i=1:length(obj.nnr)
                obj.dr(i)=(obj.rgrid(sum(obj.nnr(1:i))+1)-obj.rgrid(rid))/double(obj.nnr(i));
                rid=rid+obj.nnr(i);
            end
            
            Br = h5read(obj.fullpath,'/data/fields/Br')*obj.bnorm;
            obj.Br= reshape(Br,length(obj.zgrid),length(obj.rgrid));
            Bz = h5read(obj.fullpath,'/data/fields/Bz')*obj.bnorm;
            obj.Bz= reshape(Bz,length(obj.zgrid),length(obj.rgrid));
            try
                Atheta = h5read(obj.fullpath,'/data/fields/Athet')*obj.bnorm;
                obj.Athet= reshape(Atheta,length(obj.zgrid),length(obj.rgrid));
                [rmeshgrid,~]=meshgrid(obj.rgrid,obj.zgrid);
                obj.rAthet=(rmeshgrid.*obj.Athet)';
            catch
            end
            obj.B=sqrt(obj.Bz.^2+obj.Br.^2);
            obj.costhet=(obj.Br./obj.B)';
            obj.sinthet=(obj.Bz./obj.B)';
            clear Br Bz
            try
                obj.t0d=h5read(obj.fullpath,'/data/var0d/time');
            catch
                obj.t0d=obj.dt.*double(0:length(obj.epot)-1);
            end
            try
                for i=0:nbrst
                    grp=sprintf('/data/input.%02i/',i);
                    obj.Erxt(:,:,i+1)=reshape(h5read(obj.fullpath,[grp,'Erxt']),length(obj.zgrid),length(obj.rgrid))'*obj.enorm;
                    obj.Ezxt(:,:,i+1)=reshape(h5read(obj.fullpath,[grp,'Ezxt']),length(obj.zgrid),length(obj.rgrid))'*obj.enorm;
                    obj.potxt(:,:,i+1)=reshape(h5read(obj.fullpath,[grp,'potxt']),length(obj.zgrid),length(obj.rgrid))'*obj.phinorm;
                end
            catch
            end
            
            
            obj.femorder = h5read(obj.fullpath,'/data/input.00/femorder');
            obj.ngauss = h5read(obj.fullpath,'/data/input.00/ngauss');
            obj.plasmadim = h5read(obj.fullpath,'/data/input.00/plasmadim');
            obj.radii = h5read(obj.fullpath,'/data/input.00/radii');
            
            obj.epot = h5read(obj.fullpath,'/data/var0d/epot');
            obj.ekin = h5read(obj.fullpath,'/data/var0d/ekin');
            obj.etot = h5read(obj.fullpath,'/data/var0d/etot');
            try
                obj.etot0 = h5read(obj.fullpath,'/data/var0d/etot0');
                obj.eerr = obj.etot-obj.etot0;
            catch
                obj.eerr = obj.etot-obj.etot(2);
            end
            
            if(obj.normalized)
                obj.pot=gridquantity(obj.fullpath,'/data/fields/pot',sum(obj.nnr)+1, obj.nz+1,1);
                obj.Er=gridquantity(obj.fullpath,'/data/fields/Er',sum(obj.nnr)+1, obj.nz+1,1);
                obj.Ez=gridquantity(obj.fullpath,'/data/fields/Ez',sum(obj.nnr)+1, obj.nz+1,1);
            else
                obj.pot=gridquantity(obj.fullpath,'/data/fields/pot',sum(obj.nnr)+1, obj.nz+1,obj.phinorm);
                obj.Er=gridquantity(obj.fullpath,'/data/fields/Er',sum(obj.nnr)+1, obj.nz+1,obj.enorm);
                obj.Ez=gridquantity(obj.fullpath,'/data/fields/Ez',sum(obj.nnr)+1, obj.nz+1,obj.enorm);
            end
            
            try
                obj.t2d = h5read(obj.fullpath,'/data/fields/time');
            catch
                info=h5info(obj.fullpath,'/data/fields/partdensity');
                obj.t2d=obj.dt*(0:info.objspace.Size(2)-1)*double(obj.it1);
            end
            
            try
                info=h5info(obj.fullpath,'/data/fields/moments');
                obj.femorder = h5read(obj.fullpath,'/data/input.00/femorder');
                kr=obj.femorder(2)+1;
                obj.knotsr=augknt(obj.rgrid,kr);
                
                kz=obj.femorder(1)+1;
                obj.knotsz=augknt(obj.zgrid,kz);
                try
                    obj.CellVol= reshape(h5read(obj.fullpath,'/data/fields/volume'),length(obj.knotsz)-kz,length(obj.knotsr)-kr);
                    obj.CellVol=permute(obj.CellVol,[2,1,3])*obj.rnorm^3;
                catch
                    zvol=fnder(spmak(obj.knotsz,ones(1,length(obj.knotsz)-kz)), -1 );
                    rvol=fnder(spmak(obj.knotsr,2*pi*[obj.rgrid' 2*obj.rgrid(end)-obj.rgrid(end-1)]), -1 );
                    ZVol=diff(fnval(zvol,obj.knotsz));
                    RVol=diff(fnval(rvol,obj.knotsr));
                    obj.CellVol=RVol(3:end-1)*ZVol(3:end-1)';
                    obj.CellVol=padarray(obj.CellVol,[1,1],'replicate','post');
                end
                
                try
                    obj.geomweight = h5read(obj.fullpath,'/data/input.00/geometry/geomweight');
                    obj.geomweight= reshape(obj.geomweight,length(obj.zgrid),length(obj.rgrid),[]);
                    obj.geomweight = permute(obj.geomweight,[2,1,3]);
                catch
                    obj.geomweight=ones(length(obj.rgrid),length(obj.zgrid),3);
                end
                try
                    obj.dirichletweight = h5read(obj.fullpath,'/data/input.00/geometry/dirichletweight');
                    obj.dirichletweight= reshape(obj.dirichletweight,length(obj.zgrid),length(obj.rgrid),[]);
                    obj.dirichletweight = permute(obj.dirichletweight,[2,1,3]);
                catch
                    obj.dirichletweight=obj.geomweight;
                end
                try
                    obj.gtilde = h5read(obj.fullpath,'/data/input.00/geometry/gtilde');
                    obj.gtilde= reshape(obj.gtilde,length(obj.zgrid),length(obj.rgrid),[]);
                    obj.gtilde = permute(obj.gtilde,[2,1,3]);
                catch
                    obj.gtilde=zeros(length(obj.rgrid),length(obj.zgrid),3);
                end
                geomweight=ones(length(obj.rgrid),length(obj.zgrid));
                if(obj.normalized)
                    obj.N=splinedensity(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.CellVol, 1, geomweight, 1);
                    obj.phi=splinequantity(obj.fullpath,'/data/fields/phi', obj.knotsr, obj.knotsz, obj.femorder, 1, obj.geomweight(:,:,1), -1);
                else
                    obj.N=splinedensity(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.CellVol, abs(obj.qsim/obj.qe), geomweight, 1);
                end
                obj.fluidUR=splinevelocity(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.vnorm, geomweight, 2);
                obj.fluidUTHET=splinevelocity(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.vnorm, geomweight, 3);
                obj.fluidUZ=splinevelocity(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.vnorm, geomweight, 4);
                if(obj.normalized)
                    obj.Presstens=splinepressure(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.CellVol, obj.vnorm^2*obj.me, geomweight);
                    obj.fluidEkin=splineenergy(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.CellVol, obj.vnorm^2*obj.me*0.5, geomweight);
                else
                    obj.Presstens=splinepressure(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.CellVol, obj.vnorm^2*obj.msim, geomweight);
                    obj.fluidEkin=splineenergy(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder, obj.CellVol, obj.vnorm^2*obj.msim*0.5, geomweight);
                end
                try
                    obj.celltype=h5read(obj.fullpath,'/data/input.00/geometry/ctype')';
                    obj.linked_s=h5read(obj.fullpath,'/data/input.00/geometry/linked_s');
                    obj.bsplinetype=h5read(obj.fullpath,'/data/input.00/geometry/bsplinetype');
                    obj.bsplinetype=reshape(obj.bsplinetype,length(obj.knotsz)-kz,length(obj.knotsr)-kr);
                catch
                    obj.celltype=[];
                    obj.linked_s=[];
                end
            catch
                obj.CellVol=(obj.zgrid(2:end)-obj.zgrid(1:end-1))*((obj.rgrid(2:end).^2-obj.rgrid(1:end-1).^2)*pi)';
                obj.CellVol=obj.CellVol';
                obj.N=griddensity(obj.fullpath, '/data/fields/partdensity', sum(obj.nnr)+1, obj.nz+1, obj.CellVol, abs(obj.qsim/obj.qe), true);
                obj.fluidUR=gridquantity(obj.fullpath, '/data/fields/fluidur', sum(obj.nnr)+1, obj.nz+1, obj.vnorm, true);
                obj.fluidUTHET=gridquantity(obj.fullpath, '/data/fields/fluiduthet', sum(obj.nnr)+1, obj.nz+1, obj.vnorm, true);
                obj.fluidUZ=gridquantity(obj.fullpath, '/data/fields/fluiduz', sum(obj.nnr)+1, obj.nz+1, obj.vnorm, true);
            end
            
            % If we have a maxwellian source, read its parameters
            try
                obj.maxwellsrce.rlim=h5read(obj.fullpath, '/data/input.00/maxwellsource/rlimits');
                obj.maxwellsrce.zlim=h5read(obj.fullpath, '/data/input.00/maxwellsource/zlimits');
                obj.maxwellsrce.frequency=h5readatt(obj.fullpath, '/data/input.00/maxwellsource','frequency');
                obj.maxwellsrce.radialtype=h5readatt(obj.fullpath, '/data/input.00/maxwellsource','radialtype');
                obj.maxwellsrce.temperature=h5readatt(obj.fullpath, '/data/input.00/maxwellsource','temperature');
                obj.maxwellsrce.time_end=h5readatt(obj.fullpath, '/data/input.00/maxwellsource','time_end');
                obj.maxwellsrce.time_start=h5readatt(obj.fullpath, '/data/input.00/maxwellsource','time_start');
                obj.maxwellsrce.vth=h5readatt(obj.fullpath, '/data/input.00/maxwellsource','vth');
                obj.maxwellsrce.rate=obj.maxwellsrce.frequency*obj.weight/(pi*(diff(obj.maxwellsrce.rlim.^2))*diff(obj.maxwellsrce.zlim));
                obj.maxwellsrce.current=obj.maxwellsrce.frequency*obj.weight*obj.qe;
                obj.maxwellsrce.present=true;
            catch
                obj.maxwellsrce.present=false;
            end
            %% load neutcol parameters
            try
                obj.neutcol.neutdens=double(h5readatt(obj.fullpath, '/data/input.00/neutcol','neutdens'));
                obj.neutcol.neutpressure=double(h5readatt(obj.fullpath, '/data/input.00/neutcol','neutpressure'));
                obj.neutcol.scatter_fac=double(h5readatt(obj.fullpath, '/data/input.00/neutcol','scatter_fac'));
                obj.neutcol.Eion=double(h5readatt(obj.fullpath, '/data/input.00/neutcol','Eion'));
                obj.neutcol.E0=double(h5readatt(obj.fullpath, '/data/input.00/neutcol','E0'));
                obj.neutcol.Escale=double(h5readatt(obj.fullpath, '/data/input.00/neutcol','Escale'));
                try
                    obj.neutcol.io_cross_sec=double(h5read(obj.fullpath, '/data/input.00/neutcol/io_cross_sec'));
                    obj.neutcol.io_cross_sec(:,2)=obj.neutcol.io_cross_sec(:,2)*obj.rnorm^2;
                    obj.neutcol.io_cross_sec(:,3)=[log(obj.neutcol.io_cross_sec(2:end,2)./obj.neutcol.io_cross_sec(1:end-1,2))...
                        ./log(obj.neutcol.io_cross_sec(2:end,1)./obj.neutcol.io_cross_sec(1:end-1,1)); 0];
                    obj.neutcol.iom_cross_sec=zeros(500,3);
                    obj.neutcol.iom_cross_sec(:,1)=logspace(log10(obj.neutcol.Eion+0.001),log10(5e4),size(obj.neutcol.iom_cross_sec,1));
                    obj.neutcol.iom_cross_sec(:,2)=obj.sigmiopre(obj.neutcol.iom_cross_sec(:,1),true);
                    obj.neutcol.iom_cross_sec(:,3)=abs([log(obj.neutcol.iom_cross_sec(2:end,2)./obj.neutcol.iom_cross_sec(1:end-1,2))...
                        ./log(obj.neutcol.iom_cross_sec(2:end,1)./obj.neutcol.iom_cross_sec(1:end-1,1)); 0]);
                    
                catch
                    obj.neutcol.io_cross_sec=[];
                    obj.neutcol.iom_cross_sec=[];
                end
                try
                    obj.neutcol.ela_cross_sec=double(h5read(obj.fullpath, '/data/input.00/neutcol/ela_cross_sec'));
                    obj.neutcol.ela_cross_sec(:,2)=obj.neutcol.ela_cross_sec(:,2)*obj.rnorm^2;
                    obj.neutcol.ela_cross_sec(:,3)=[log(obj.neutcol.ela_cross_sec(2:end,2)./obj.neutcol.ela_cross_sec(1:end-1,2))...
                        ./log(obj.neutcol.ela_cross_sec(2:end,1)./obj.neutcol.ela_cross_sec(1:end-1,1)); 0];
                catch
                    obj.neutcol.ela_cross_sec=[];
                end
                obj.neutcol.present=true;
            catch
                obj.neutcol.present=false;
            end
            
            %% load spline boundaries
            try
                obj.spl_bound.nbsplines=h5readatt(obj.fullpath, '/data/input.00/geometry_spl','nbsplines');
                for i=1:obj.spl_bound.nbsplines
                    splgroup=sprintf('/data/input.00/geometry_spl/%02d',i);
                    obj.spl_bound.boundary(i).knots=h5read(obj.fullpath,sprintf('%s/knots',splgroup));
                    obj.spl_bound.boundary(i).Dval=h5readatt(obj.fullpath,splgroup,'Dirichlet_val');
                    obj.spl_bound.boundary(i).coefs=reshape(h5read(obj.fullpath,sprintf('%s/pos',splgroup)),2,[])';
                    obj.spl_bound.boundary(i).order=h5readatt(obj.fullpath,splgroup,'order');
                    obj.spl_bound.boundary(i).kind=h5readatt(obj.fullpath,splgroup,'kind');
                    obj.spl_bound.boundary(i).fun=spmak(obj.spl_bound.boundary(i).knots,obj.spl_bound.boundary(i).coefs');
                end
            catch
                obj.spl_bound.nbsplines=0;
            end
            
             %% load non ideal power supply parameters
            try
                obj.psupply.targetbias=h5readatt(obj.fullpath, '/data/input.00/psupply','targetbias');
                obj.psupply.expdens=h5readatt(obj.fullpath, '/data/input.00/psupply','expdens');
                obj.psupply.PSresistor=h5readatt(obj.fullpath, '/data/input.00/psupply','PSresistor');
                obj.psupply.geomcapacitor=h5readatt(obj.fullpath, '/data/input.00/psupply','geomcapacitor');
                obj.psupply.nbhdt=h5readatt(obj.fullpath, '/data/input.00/psupply','nbhdt');
                obj.psupply.biases=h5read(obj.fullpath, '/data/var0d/biases');
                obj.psupply.current=h5read(obj.fullpath, '/data/var0d/current');
                obj.psupply.tau=obj.psupply.PSresistor*obj.psupply.geomcapacitor*obj.psupply.expdens/obj.neutcol.neutdens;
                obj.psupply.active=true;
                obj.psupply.bdpos=h5read(obj.fullpath, '/data/input.00/psupply/bdpos');
            catch
                obj.psupply.active=false;
            end
            
            % Read the main particles parameters
            if(readparts)
                
                if(obj.normalized)
                    obj.R = h5partsquantity(obj.fullpath,'/data/part','R');
                    obj.Z = h5partsquantity(obj.fullpath,'/data/part','Z');
                else
                    obj.R = h5partsquantity(obj.fullpath,'/data/part','R',obj.rnorm);
                    obj.Z = h5partsquantity(obj.fullpath,'/data/part','Z',obj.rnorm);
                end
                try
                    obj.THET = h5partsquantity(obj.fullpath,'/data/part','THET');
                catch
                    clear obj.THET
                end
                try
                    obj.Rindex=h5partsquantity(obj.fullpath,'/data/part','Rindex');
                    obj.Zindex=h5partsquantity(obj.fullpath,'/data/part','Zindex');
                catch
                    clear obj.Rindex obj.Zindex
                end
                vscale=obj.vnorm;
                
                obj.VR = h5partsquantity(obj.fullpath,'/data/part','UR',vscale);
                obj.VZ = h5partsquantity(obj.fullpath,'/data/part','UZ',vscale);
                obj.VTHET= h5partsquantity(obj.fullpath,'/data/part','UTHET',vscale);
                
                if(obj.normalized)
                    obj.partepot = h5partsquantity(obj.fullpath,'/data/part','pot',sign(obj.qsim)*obj.qe);
                else
                    obj.partepot = h5partsquantity(obj.fullpath,'/data/part','pot',sign(obj.qsim)*obj.qe*obj.phinorm);
                end
                
                try
                    obj.partindex = h5partsquantity(obj.fullpath,'/data/part/','partindex');
                catch
                end
                
                if(obj.nbspecies >1)
                    obj.species=h5parts.empty(obj.nbspecies,0);
                    for i=2:obj.nbspecies
                        obj.species(i-1)=h5parts(obj.fullpath,sprintf('/data/part/%2d',i),obj);
                    end
                end
            end
            
            try
                obj.tpart = h5read(obj.fullpath,'/data/part/time');
                obj.nbparts = h5read(obj.fullpath,'/data/part/Nparts');
            catch
                obj.nbparts=obj.npart;
                obj.tpart=obj.dt*(0:size(obj.R,2)-1)*double(obj.it2);
            end
            
            
            if(obj.nbcelldiag > 0)
                obj.celldiag=h5parts.empty;
                j=0;
                for i=1:obj.nbcelldiag
                    nbparts=h5read(obj.fullpath,sprintf('%s/Nparts',sprintf('/data/celldiag/%02d',i)));
                    if (sum(nbparts)>0)
                        j=j+1;
                        obj.celldiag(j)=h5parts(obj.fullpath,sprintf('/data/celldiag/%02d',i),obj);
                        obj.celldiag(j).rindex=double(h5readatt(obj.fullpath, sprintf('/data/celldiag/%02d',i),'rindex'))+(1:2);
                        obj.celldiag(j).zindex=double(h5readatt(obj.fullpath, sprintf('/data/celldiag/%02d',i),'zindex'))+(1:2);
                    end
                    
                end
            end
        end
        %------------------------------------------
        %  Functions for accesing secondary simulation quantities
        function Atheta=Atheta(obj,R,Z)
            %% returns the magnetic vector potential at position R,Z interpolated from stored Athet in h5 file
            %             halflz=(obj.zgrid(end)+obj.zgrid(1))/2;
            %             Atheta=0.5*obj.B0*(R-obj.width/pi*(obj.Rcurv-1)/(obj.Rcurv+1)...
            %                 .*besseli(1,2*pi*R/obj.width).*cos(2*pi*(Z-halflz)/obj.width));
            Atheta=interp2(obj.rgrid,obj.zgrid,obj.Athet,R,Z);
        end
        
        function quantity=H(obj,indices)
            %% computes the total energy for the main specie
            % for the particle with index indices{1} at timepart step indices{2}
            % which is time obj.timepart(indices{2})
            if strcmp(indices{1},':')
                p=1:obj.VR.nparts;% if nothing is defined we load all particles
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart); %if nothing is defined all time steps are considered
            else
                t=indices{2};
            end
            % if track is true we look at specific particles with their
            % index and follow them in time
            % if it is false we just care about the distribution function
            % and specific particles can have different positions in the
            % resulting array for each timestep
            if size(indices,1)>2
                track=indices{3};
            else
                track=false;
            end
            quantity=0.5*obj.me*(obj.VR(p,t,track).^2+obj.VTHET(p,t,track).^2+obj.VZ(p,t,track).^2)+obj.partepot(p,t,track);
        end
        
        function quantity=P(obj,indices)
            %P computes the canonical angular momentum for the main specie
            % for the particle with index indices{1} at timepart step indices{2}
            % which is time obj.timepart(indices{2})
            if strcmp(indices{1},':')
                p=1:obj.R.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            % if track is true we look at specific particles with their
            % index and follow them in time
            % if it is false we just care about the distribution function
            % and specific particles can have different positions in the
            % resulting array for each timestep
            if size(indices,1)>2
                track=indices{3};
            else
                track=false;
            end
            quantity=obj.R(p,t,track).*(obj.VTHET(p,t,track)*obj.me+sign(obj.qsim)*obj.qe*obj.Atheta(obj.R(p,t,track),obj.Z(p,t,track)));
        end
        
        function quantity=Vpar(obj,varargin)
            %Vpar Computes the parallel velocity for the main specie
            % for the particle with index indices{1} at timepart step indices{2}
            % which is time obj.timepart(indices{2})
            if(~iscell(varargin))
                indices=mat2cell(varargin);
            else
                indices=varargin;
            end
            if strcmp(indices{1},':')
                p=1:obj.R.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            % if track is true we look at specific particles with their
            % index and follow them in time
            % if it is false we just care about the distribution function
            % and specific particles can have different positions in the
            % resulting array for each timestep
            if size(indices,1)>2
                track=indices{3};
            else
                track=false;
            end
            Zp=obj.Z(p,t,track);% get the particle axial positon
            Rp=obj.R(p,t,track);% get the particle radial position
            
            % interpolate the magnetic field at the particle position
            Bzp=interp2(obj.zgrid,obj.rgrid,obj.Bz',Zp,Rp,'makima');
            Brp=interp2(obj.zgrid,obj.rgrid,obj.Br',Zp,Rp,'makima');
            Bp=sqrt(Bzp.^2+Brp.^2);
            % calculate the projection angle of the radial and axial
            % directions on the magnetic field line
            Costhet=Bzp./Bp;
            Sinthet=Brp./Bp;
            % calculate the actuale parallel velocity
            quantity=obj.VR(p,t,track).*Sinthet+obj.VZ(p,t,track).*Costhet;
        end
        
        function quantity=Vperp(obj,varargin)
            %Vperp Computes the perpendicular velocity in the guidind center reference frame,
            % for the main specie particle indices{1} at time indices{2}
            
            if(~iscell(varargin))
                indices=mat2cell(varargin);
            else
                indices=varargin;
            end
            
            if strcmp(indices{1},':')
                p=1:obj.R.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            % if track is true we look at specific particles with their
            % index and follow them in time
            % if it is false we just care about the distribution function
            % and specific particles can have different positions in the
            % resulting array for each timestep
            if size(indices,2)>2
                track=indices{3};
            else
                track=false;
            end
            
            % if gcs is true, gives the perpendicular velocity in the
            % guiding center system by substracting the EXB azimuthal
            % velocity
            % else gives the total perpendicular velocity
            if size(indices,2)>3
                gcs=indices{4};
            else
                gcs=false;
            end
            % get the particle position
            Zp=obj.Z(p,t,track);
            Rp=obj.R(p,t,track);
            % interpolate the magnetic field at the particle position
            Bzp=interp2(obj.zgrid,obj.rgrid,obj.Bz',Zp,Rp,'makima');
            Brp=interp2(obj.zgrid,obj.rgrid,obj.Br',Zp,Rp,'makima');
            Bp=sqrt(Bzp.^2+Brp.^2);
            % calculate the projecting angles
            Costhet=Bzp./Bp;
            Sinthet=Brp./Bp;
            Vdrift=zeros(size(Zp));
            
            if gcs
                % for each particle and each timestep
                % calculate the azimuthal ExB drift velocity
                for j=1:length(t)
                    [~, tfield]=min(abs(obj.t2d-obj.tpart(t(j))));
                    timeEr=obj.Er(:,:,tfield);
                    timeEz=obj.Ez(:,:,tfield);
                    %posindE=sub2ind(size(timeEr),Rind(:,j),Zind(:,j));
                    timeErp=interp2(obj.zgrid,obj.rgrid,timeEr,Zp(:,j),Rp(:,j));
                    timeEzp=interp2(obj.zgrid,obj.rgrid,timeEz,Zp(:,j),Rp(:,j));
                    Vdrift(:,j)=(timeEzp.*Brp(:,j)-timeErp.*Bzp(:,j))./Bp(:,j).^2;
                end
            end
            % calculate the perpendicular velocity
            quantity=sqrt((obj.VTHET(p,t,track)-Vdrift).^2+(obj.VR(p,t,track).*Costhet-obj.VZ(p,t,track).*Sinthet).^2);
        end
        
        function quantity=cyclphase(obj,varargin)
            %cyclphase Computes the cyclotronic phase for the main specie
            % for particles with indices{1} at time indices{2}
            
            if(~iscell(varargin))
                indices=mat2cell(varargin);
            else
                indices=varargin;
            end
            
            if strcmp(indices{1},':')
                p=1:obj.R.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            % if track is true we look at specific particles with their
            % index and follow them in time
            % if it is false we just care about the distribution function
            % and specific particles can have different positions in the
            % resulting array for each timestep
            if size(indices,2)>2
                track=indices{3};
            else
                track=false;
            end
            Zp=obj.Z(p,t,track);
            Rp=obj.R(p,t,track);
            %             [~, zind(1)]=min(abs(obj.zgrid-0.005262));
            %             [~, zind(2)]=min(abs(obj.zgrid-0.006637));
            %             [~, rind(1)]=min(abs(obj.rgrid-0.0784));
            %             [~, rind(2)]=min(abs(obj.rgrid-0.07861));
            %             indices=Zp<obj.zgrid(zind(2)) & Zp>=obj.zgrid(zind(1)) &...
            %                 Rp<obj.rgrid(rind(2)) & Rp>=obj.rgrid(rind(1));
            %Zp=Zp(indices);
            %Rp=Rp(indices);
            %p=indices;
            
            Bzp=interp2(obj.zgrid,obj.rgrid,obj.Bz',Zp,Rp,'makima');
            Brp=interp2(obj.zgrid,obj.rgrid,obj.Br',Zp,Rp,'makima');
            Bp=sqrt(Bzp.^2+Brp.^2);
            Costhet=Bzp./Bp;
            Sinthet=Brp./Bp;
            
            % compute the projection  of the perpendicular velocity in the
            % radial direction
            vr=(obj.VR(p,t,track).*Costhet-obj.VZ(p,t,track).*Sinthet);
            % Get the perpendicular velocity
            vperp=obj.Vperp(p,t,track,true);
            vr=vr(indices);
            vperp=vperp(indices);
            cospsi=vr./vperp;
            quantity=acos(cospsi);
        end
        
        function p=borderpoints(obj,subdiv)
            %borderpoints Return a cell array containing the curves
            %defining the boundary of the domain
            % for each boundary p(1,:) and p(2,:) give axial and radial position
            % for each boundary p(3,:) and p(4,:) give axial and radial normals
            
            %gw= contourc(obj.zgrid,obj.rgrid,obj.geomweight(:,:,1),[0 0])
            p=cell(0,0);
            if nargin<2
                subdiv=1;
            end
            
            ndiv=sum(subdiv);
            
            %outer cylinder
            if any(obj.geomweight(end,:,1)>=0)
                idp=ceil(length(obj.zgrid)/ndiv);
                imin=1;
                for j=1:length(subdiv)
                    imax=min(imin+subdiv(j)*idp-1,length(obj.zgrid));
                    p{end+1}=[obj.zgrid(imin:imax)';obj.rgrid(end)*ones(imax-imin+1,1)';
                        zeros(imax-imin+1,1)';ones(imax-imin+1,1)'];
                    imin=imax;
                end
                
            end
            
            %inner cylinder
            if any(obj.geomweight(1,:,1)>=0)
                idp=ceil(length(obj.zgrid)/ndiv);
                imin=1;
                for j=1:length(subdiv)
                    imax=min(imin+subdiv(j)*idp-1,length(obj.zgrid));
                    p{end+1}=[obj.zgrid(imin:imax)';obj.rgrid(1)*ones(imax-imin+1,1)';
                        zeros(imax-imin+1,1)';-ones(imax-imin+1,1)'];
                    imin=imax;
                end
            end
            
            if obj.walltype==2
                % We have an elliptic insert that we want to isolate
                gw=obj.ellipseborder;
                zpos=obj.zgrid(obj.zgrid<(min(gw(1,:))) | obj.zgrid>(max(gw(1,:))));
                p{2}=[zpos,obj.rgrid(end)*ones(size(zpos))]';
                p{1}=[obj.zgrid';obj.rgrid(1)*ones(size(obj.zgrid))'];
                gw=obj.ellipseborder;
                p{3}=gw;
            elseif obj.walltype~=0
                % extract all the walls
                gw=contourc(obj.zgrid,obj.rgrid,obj.geomweight(:,:,1),[0 0]);
                [x,y,~]=C2xyz(gw);
                for i=1:length(x)
                    %subdiv=[4,1,2];
                    ndiv=sum(subdiv);
                    idp=ceil(length(x{i})/ndiv);
                    imin=1;
                    for j=1:length(subdiv)
                        imax=min(imin+subdiv(j)*idp-1,length(x{i}));
                        p{end+1}=[x{i}(imin:imax);y{i}(imin:imax)];
                        imin=imax;
                    end
                end
            end
            %             figure
            %             for i=1:length(p)
            %                 plot(p{i}(1,:),p{i}(2,:))
            %                 hold on
            %             end
        end
        
        function p=ellipseborder(obj)
            %ellipseborder returns the boundary points defining the
            %elliptic insert
            z=linspace(-0.998,0.998,1000)*obj.z_r;
            p=zeros(4,length(z));
            for i=1:length(z)
                p(1,i)=z(i)+obj.z_0;
                p(2,i)=obj.r_0-obj.r_r*sqrt(1-(z(i)/obj.z_r)^2);
                p(3,i)=2/(obj.z_r^2)*(z(i));
                p(4,i)=2/(obj.r_r^2)*(p(2,i)-obj.r_0);
            end
            norm=sqrt(p(3,:).^2+p(4,:).^2);
            p(3,:)=double(obj.interior)*p(3,:)./norm;
            p(4,:)=double(obj.interior)*p(4,:)./norm;
        end
        
        function charge=totcharge(obj,fieldstep)
            % Integrates the density profile over the full volume to obtain
            % the total number of electrons in the volume
            n=splinedensity(obj.fullpath, '/data/fields/moments', obj.knotsr, obj.knotsz, obj.femorder,ones(size(obj.CellVol)), 1, 1);
            charge=sum(sum(n(:,:,fieldstep)));
        end
        
        function Gamma=Axialflux(obj,timestep,zpos)
            % Computes the axial particle flux n*Uz at timestep timestep and axial position zpos
            Gamma=obj.fluidUZ(:,zpos,timestep).*obj.N(:,zpos,timestep);
        end
        
        function Gamma=Metallicflux(obj,timestep,subdiv)
            % Computes the particle flux at time obj.t2d(timestep) on the
            % metallic boundaries
            
            if nargin<3
                subdiv=1;
            end
            % We find the borderpoints
            p=obj.borderpoints(subdiv);
            gamma=cell(size(p));
            Nr=cell(size(p));
            Nz=cell(size(p));
            for i=1:length(p)
                bp=p{i};
                if size(bp,1)==2
                    % We get the normals at these positions and normalise them
                    Nr{i}=-interp2(obj.zgrid,obj.rgrid,obj.geomweight(:,:,3),bp(1,:),bp(2,:));
                    Nz{i}=-interp2(obj.zgrid,obj.rgrid,obj.geomweight(:,:,2),bp(1,:),bp(2,:));
                    norm=sqrt(Nr{i}.^2+Nz{i}.^2);
                    Nr{i}=Nr{i}./norm;
                    Nz{i}=Nz{i}./norm;
                else
                    Nr{i}=bp(4,:);
                    Nz{i}=bp(3,:);
                end
                
                gamma{i}=zeros(size(bp,2),length(timestep));
            end
            [z,r]=ndgrid(obj.zgrid,obj.rgrid);
            N=obj.N(:,:,timestep(1));
            n=griddedInterpolant(z,r,N');
            uz=griddedInterpolant(z,r,obj.fluidUZ(:,:,timestep(1))');
            ur=griddedInterpolant(z,r,obj.fluidUR(:,:,timestep(1))');
            % we get the density and fluid velocities at the desired time
            % steps and interpolate them at the boundary position
            for j=1:length(timestep)
                n.Values=obj.N(:,:,timestep(j))';
                uz.Values=obj.fluidUZ(:,:,timestep(j))';
                ur.Values=obj.fluidUR(:,:,timestep(j))';
                for i=1:length(p)
                    bp=p{i};
                    gamma{i}(:,j)=n(bp(1:2,:)').*(ur(bp(1:2,:)').*Nr{i}'+uz(bp(1:2,:)').*Nz{i}');
                end
            end
            % return the boundary position p and the corresponding flux
            % gamma
            Gamma.p=p;
            Gamma.gamma=gamma;
        end
        
        function [I, pos]=OutCurrents(obj,timestep, subdiv)
            % Computes the Outgoing currens at the simulation axial boundaries at timestep timestep
            % This is simply the surface integral of the axial flux
            if nargin<3
                subdiv=1;
            end
            flux=obj.Axialflux(timestep,[1 obj.nz+1]);
            Iz=squeeze(trapz(obj.rgrid,flux.*obj.rgrid)*2*pi*obj.qsim/obj.weight);
            Iz(1,:)=-Iz(1,:);
            gamm=obj.Metallicflux(timestep, subdiv);
            Im=zeros(length(gamm.p),length(timestep));
            pos=cell(size(gamm.p));
            for i=1:length(gamm.p)
                p=gamm.p{i};
                pos{i}=p;
                flux=gamm.gamma{i};
                for j=1:length(timestep)
                    Im(i,j)=pi/2*sum((p(2,1:end-1)+p(2,2:end)).*(flux(2:end,j)+flux(1:end-1,j))'...
                        .*sqrt((p(1,2:end)-p(1,1:end-1)).^2+(p(2,2:end)-p(2,1:end-1)).^2));
                end
            end
            I=-cat(1,Iz,Im*obj.qsim/obj.weight);
        end
        
        function [pot] = PotentialWell(obj,fieldstep)
            %PotentialWell Computes the potential well at the given timestep on the FEM grid points
            % interpolates the model data on rgrid and zgrid
            model=obj.potentialwellmodel(fieldstep);
            z=model.z;
            r=model.r;
            modpot=model.pot;
            [Zmesh,Rmesh]=meshgrid(obj.zgrid,obj.rgrid);
            pot=zeros(length(obj.zgrid),length(obj.rgrid),length(fieldstep));
            for i=1:length(fieldstep)
                pot(:,:,i)=griddata(z,r,modpot(:,i),Zmesh,Rmesh)';
            end
        end
        
        function Epar = Epar(obj,fieldstep)
            % Computes the electric field component parallel to the magnetic field line
            Epar=obj.Er(:,:,fieldstep).*(obj.Br./obj.B)' + (obj.Bz./obj.B)'.*obj.Ez(:,:,fieldstep);
        end
        
        function Eperp = Eperp(obj,fieldstep)
            % Computes the electric field component perpendicular to the magnetic field line
            Eperp=obj.Er(:,:,fieldstep).*(obj.Bz./obj.B)' - (obj.Br./obj.B)'.*obj.Ez(:,:,fieldstep);
        end
        
        function Ekin = Ekin(obj,varargin)
            %Ekin Computes the classical kinetic energy of particles indices{1} at
            % time obj.tpart(indices{2}) in Joules
            if(~iscell(varargin))
                indices=mat2cell(varargin);
            else
                indices=varargin;
            end
            if strcmp(indices{1},':')
                p=1:obj.R.nparts;
            else
                p=indices{1};
            end
            if strcmp(indices{2},':')
                t=1:length(obj.tpart);
            else
                t=indices{2};
            end
            % if track is true we look at specific particles with their
            % index and follow them in time
            % if it is false we just care about the distribution function
            % and specific particles can have different positions in the
            % resulting array for each timestep
            if size(indices,1)>2
                track=indices{3};
            else
                track=false;
            end
            Vr=obj.VR(p,t,track);
            Vthet= obj.VTHET(p,t,track);
            Vz=obj.VZ(p,t,track);
            Ekin=0.5*obj.msim/obj.weight*(Vr.^2+Vthet.^2+Vz.^2);
        end
        
        function sig=sigio(obj,E,init)
            %sigio returns the total ionisation cross-section in m^2
            % at energy E[eV]
            % init is only used during the loading of the h5 file
            if nargin <3
                init=false;
            end
            sig=zeros(size(E));
            
            if(~init &&( ~obj.neutcol.present || isempty(obj.neutcol.io_cross_sec)))
                sig=zeros(size(E));
                return
            end
            for i=1:length(E(:))
                if(E(i)>obj.neutcol.Eion)
                    sig(ind2sub(size(E),i))=obj.fit_cross_sec(E(ind2sub(size(E),i)),obj.neutcol.io_cross_sec);
                end
            end
        end
        
        function sig=sigmio(obj,E)
            %sigmio returns the total ionisation cross-section for momentum exchange for the incoming electron in m^2
            % at energy E[eV]
            sig=zeros(size(E));
            if(~obj.neutcol.present || isempty(obj.neutcol.iom_cross_sec))
                return
            end
            for i=1:length(E(:))
                if(E(i)>obj.neutcol.Eion)
                    sig(ind2sub(size(E),i))=obj.fit_cross_sec(E(ind2sub(size(E),i)),obj.neutcol.iom_cross_sec);
                end
            end
        end
        
        function sigm=sigmela(obj,E)
            %sigmela returns the elastic collision cross-section for momentum exchange for the incoming electron in m^2
            % at energy E[eV]
            sigm=zeros(size(E));
            if(~obj.neutcol.present || isempty(obj.neutcol.ela_cross_sec))
                return
            end
            for i=1:length(E(:))
                sigm(ind2sub(size(E),i))=obj.fit_cross_sec(E(ind2sub(size(E),i)),obj.neutcol.ela_cross_sec);
            end
        end
        
        function sig=sigela(obj,E)
            %sigmela returns the elastic collision cross-section for the incoming electron in m^2
            % at energy E[eV]
            % if used this will give the frequency of elastic collisions
            E0=obj.neutcol.E0;
            chi=E./(0.25*E0+E);
            sig=(2*chi.^2)./((1-chi).*((1+chi).*log((1+chi)./(1-chi))-2*chi)).*obj.sigmela(E);
        end
        
        function [Forces, Density]=Forcespline(obj,it,fdens,getmean)
            %Forcesplie calculates the fluid force terms in each direction
            %at time obj.t2d(it)
            %   if fdens return the force density in N/m^3 othewise give
            %   the force in N
            %   if getmean return only the time averaged quanties over
            %   time samples[it(1)...it(end]
            
            if strcmp(it,':')
                it=floor(0.95*size(obj.t2d)):size(obj.t2d)-1;
            end
            if nargin<3
                fdens=true;
            end
            if nargin <4
                getmean=false;
            end
            
            % To be able to calculate the centered finite difference in
            % time, we remove the first and last time indices
            it(it<2)=[];
            it(it>length(obj.t2d)-1)=[];
            
            m_e=obj.msim/obj.weight;
            q_e=obj.qsim/obj.weight;
            n=obj.N(:,:,it);
            [r,~]=meshgrid(obj.rgrid,obj.zgrid);
            Rinv=1./r';
            Rinv(isinf(Rinv))=0;
            % get inverse of density to get the force in N
            Density.N=n;
            invn=1./n;
            invn(isnan(invn) | isinf(invn))=0;
            
            
            % Calculate electric forces
            Eforcer=q_e*obj.Er(:,:,it);
            Eforcez=q_e*obj.Ez(:,:,it);
            
            
            
            Dragforcer=zeros(size(n,1),size(n,2),size(n,3));
            Dragforcethet=zeros(size(n,1),size(n,2),size(n,3));
            Dragforcez=zeros(size(n,1),size(n,2),size(n,3));
            
            time=obj.t2d(it);
            Forces.it=it;
            Forces.time=time;
            
            if getmean
                if ~fdens
                    n=ones(size(n));
                end
                % Electric forces
                Forces.Eforcer=mean(n.*q_e.*obj.Er(:,:,it),3);
                Forces.Eforcez=mean(n.*q_e.*obj.Ez(:,:,it),3);
                
                % Magnetic forces
                Forces.Bforcer=mean(q_e.*obj.fluidUTHET(:,:,it).*obj.Bz'.*n,3);
                Forces.Bforcethet=mean(q_e.*(obj.fluidUZ(:,:,it).*obj.Br'-obj.fluidUR(:,:,it).*obj.Bz').*n,3);
                Forces.Bforcez=mean(-q_e.*obj.fluidUTHET(:,:,it).*obj.Br'.*n,3);
                
                % Inertial forces
                Forces.inertforcer=mean(-m_e.*n.*(-obj.fluidUTHET(:,:,it).^2.*Rinv...
                    +obj.fluidUR(:,:,it).*obj.fluidUR.der(:,:,it,[1 0])...
                    +obj.fluidUZ(:,:,it).*obj.fluidUR.der(:,:,it,[0 1])),3);
                
                Forces.inertforcethet=mean(-m_e*n.*(obj.fluidUR(:,:,it).*obj.fluidUTHET(:,:,it).*Rinv...
                    +obj.fluidUR(:,:,it).*obj.fluidUTHET.der(:,:,it,[1 0])...
                    +obj.fluidUZ(:,:,it).*obj.fluidUTHET.der(:,:,it,[0 1])),3);
                
                Forces.inertforcez=mean(-m_e*n.*(obj.fluidUR(:,:,it).*obj.fluidUZ.der(:,:,it,[1 0])...
                    +obj.fluidUZ(:,:,it).*obj.fluidUZ.der(:,:,it,[0 1])),3);
                
                % Pressure forces
                Forces.Pressforcer=mean(-n.*( squeeze(obj.Presstens.der(1,:,:,it,[1 0]))...
                    + squeeze(obj.Presstens(1,:,:,it) - obj.Presstens(4,:,:,it)).*Rinv...
                    + squeeze(obj.Presstens.der(3,:,:,it,[0 1])))...
                    .*invn,3);
                
                Forces.Pressforcethet=mean(-n.*( squeeze(obj.Presstens.der(2,:,:,it,[1 0]))...
                    + squeeze(obj.Presstens.der(5,:,:,it,[0 1])) ...
                    + 2*squeeze(obj.Presstens(2,:,:,it)).*Rinv   ...
                    ).*invn,3);
                
                Forces.Pressforcez=mean(-n.*( squeeze(obj.Presstens.der(3,:,:,it,[1 0]))...
                    + squeeze(obj.Presstens(3,:,:,it)).*Rinv...
                    + squeeze(obj.Presstens.der(6,:,:,it,[0 1])) )...
                    .*invn,3);
                
                % ellastic coll drag forces
                if( obj.neutcol.present)
                    Ek=squeeze(obj.fluidEkin(1,:,:,it)+obj.fluidEkin(2,:,:,it)+obj.fluidEkin(3,:,:,it));
                    sigm=obj.sigmela(Ek/obj.qe)+obj.sigio(Ek/obj.qe)+obj.sigmio(Ek/obj.qe);
                    dragfreq=obj.neutcol.neutdens.*sigm.*sqrt(2*obj.weight/obj.msim*Ek);
                    Forces.Dragforcer=mean(-m_e*n.*dragfreq.*obj.fluidUR(:,:,it),3);
                    Forces.Dragforcethet=mean(-m_e*n.*dragfreq.*obj.fluidUTHET(:,:,it),3);
                    Forces.Dragforcez=mean(-m_e*n.*dragfreq.*obj.fluidUZ(:,:,it),3);
                else
                    Forces.Dragforcer=0;
                    Forces.Dragforcethet=0;
                    Forces.Dragforcez=0;
                end
                
                % effective drag frequency due to the maxwellian source
                if( obj.maxwellsrce.present)
                    dragfreqsrc=obj.maxwellsrce.frequency*obj.weight/(pi*diff(obj.maxwellsrce.zlim)*(obj.maxwellsrce.rlim(2)^2-obj.maxwellsrce.rlim(1)^2)).*invn;
                    dragfreqsrc(isinf(dragfreqsrc))=0;
                    Forces.Dragforcer=Forces.Dragforcer+mean(-n.*m_e.*dragfreqsrc.*obj.fluidUR(:,:,it),3);
                    Forces.Dragforcethet=Forces.Dragforcethet+mean(-n.*m_e.*dragfreqsrc.*obj.fluidUTHET(:,:,it),3);
                    Forces.Dragforcez=Forces.Dragforcez+mean(-n.*m_e.*dragfreqsrc.*obj.fluidUZ(:,:,it),3);
                end
                
                % Time derivative for fluid accelleration
                cdt=(obj.t2d(it+1)-obj.t2d(it-1));
                cdt=reshape(cdt,1,1,[]);
                Forces.durdt=mean(m_e*(obj.fluidUR(:,:,it+1)-obj.fluidUR(:,:,it-1))./cdt,3);
                Forces.duthetdt=mean(m_e*(obj.fluidUTHET(:,:,it+1)-obj.fluidUTHET(:,:,it-1))./cdt,3);
                Forces.duzdt=mean(m_e*(obj.fluidUZ(:,:,it+1)-obj.fluidUZ(:,:,it-1))./cdt,3);
                
            else
                % Allocate memory
                Bforcer=zeros(size(n,1),size(n,2),size(n,3));
                Bforcez=zeros(size(n,1),size(n,2),size(n,3));
                Bforcethet=zeros(size(n,1),size(n,2),size(n,3));
                inertforcer=zeros(size(n,1),size(n,2),size(n,3));
                inertforcez=zeros(size(n,1),size(n,2),size(n,3));
                inertforcethet=zeros(size(n,1),size(n,2),size(n,3));
                Pressforcer=zeros(size(n,1),size(n,2),size(n,3));
                Pressforcethet=zeros(size(n,1),size(n,2),size(n,3));
                Pressforcez=zeros(size(n,1),size(n,2),size(n,3));
                durdt=zeros(size(n,1),size(n,2),size(n,3));
                duthetdt=zeros(size(n,1),size(n,2),size(n,3));
                duzdt=zeros(size(n,1),size(n,2),size(n,3));
                fluiduThet=obj.fluidUTHET(:,:,it);
                Density.fluiduThet=fluiduThet;
                for j=1:size(n,3)
                    % Magnetic forces
                    Bforcer(:,:,j)=q_e.*fluiduThet(:,:,j).*obj.Bz';
                    Bforcethet(:,:,j)=q_e.*(obj.fluidUZ(:,:,it(j)).*obj.Br'-obj.fluidUR(:,:,it(j)).*obj.Bz');
                    Bforcez(:,:,j)=-q_e.*fluiduThet(:,:,j).*obj.Br';
                    
                    % Inertial forces
                    inertforcer(:,:,j)=-m_e.*(-fluiduThet(:,:,j).^2.*Rinv...
                        +obj.fluidUR(:,:,it(j)).*obj.fluidUR.der(:,:,it(j),[1 0])...
                        +obj.fluidUZ(:,:,it(j)).*obj.fluidUR.der(:,:,it(j),[0 1]));
                    
                    inert1=obj.fluidUR(:,:,it(j)).*fluiduThet(:,:,j).*Rinv;
                    inert2=obj.fluidUR(:,:,it(j)).*obj.fluidUTHET.der(:,:,it(j),[1 0]);
                    inert3=obj.fluidUZ(:,:,it(j)).*obj.fluidUTHET.der(:,:,it(j),[0 1]);
                    
                    inertforcethet(:,:,j)=-m_e.*(inert1...
                        +inert2...
                        +inert3);
                    
                    inertforcez(:,:,j)=-m_e.*(obj.fluidUR(:,:,it(j)).*obj.fluidUZ.der(:,:,it(j),[1 0])...
                        +obj.fluidUZ(:,:,it(j)).*obj.fluidUZ.der(:,:,it(j),[0 1]));
                    
                    % Pressure forces
                    Pr1=squeeze(obj.Presstens.der(1,:,:,it(j),[1 0]));
                    
                    Pr2=squeeze(obj.Presstens(1,:,:,it(j)) - obj.Presstens(4,:,:,it(j))).*Rinv;
                    
                    Pr3=squeeze(obj.Presstens.der(3,:,:,it(j),[0 1]));
                    
                    Pressforcer(:,:,j)=-( Pr1...
                        + Pr2...
                        + Pr3 )...
                        .*invn(:,:,j);
                    
                    Pthet1=squeeze(obj.Presstens.der(2,:,:,it(j),[1 0]));
                    Pthet2=squeeze(obj.Presstens.der(5,:,:,it(j),[0 1]));
                    Pthet3=2*squeeze(obj.Presstens(2,:,:,it(j))).*Rinv;
                    
                    Pressforcethet(:,:,j)=-( Pthet1...
                        + Pthet2 ...
                        + Pthet3   ...
                        ).*invn(:,:,j);
                    
                    Pz1=squeeze(obj.Presstens.der(3,:,:,it(j),[1 0]));
                    Pz2=squeeze(obj.Presstens(3,:,:,it(j))).*Rinv;
                    Pz3=squeeze(obj.Presstens.der(6,:,:,it(j),[0 1]));
                    Pressforcez(:,:,j)=-( Pz1...
                        + Pz2...
                        + Pz3 )...
                        .*invn(:,:,j);
                    
                    % ellastic coll drag forces
                    if( obj.neutcol.present)
                        Ek=squeeze(obj.fluidEkin(1,:,:,it(j))+obj.fluidEkin(2,:,:,it(j))+obj.fluidEkin(3,:,:,it(j)));
                        sigm=obj.sigmela(Ek/obj.qe)+obj.sigio(Ek/obj.qe)+obj.sigmio(Ek/obj.qe);
                        dragfreq=obj.neutcol.neutdens.*sigm.*sqrt(2*obj.weight/obj.msim*Ek);
                        Dragforcer(:,:,j)=-m_e*dragfreq.*obj.fluidUR(:,:,it(j));
                        Dragforcethet(:,:,j)=-m_e*dragfreq.*obj.fluidUTHET(:,:,it(j));
                        Dragforcez(:,:,j)=-m_e*dragfreq.*obj.fluidUZ(:,:,it(j));
                    end
                    
                    % effective drag frequency due to the maxwellian source
                    if( obj.maxwellsrce.present)
                        dragfreqsrc=obj.maxwellsrce.frequency*obj.weight/(pi*diff(obj.maxwellsrce.zlim)*(obj.maxwellsrce.rlim(2)^2-obj.maxwellsrce.rlim(1)^2))*invn(:,:,j);
                        dragfreqsrc(isinf(dragfreqsrc))=0;
                        Dragforcer(:,:,j)=Dragforcer(:,:,j)+-m_e*dragfreqsrc.*obj.fluidUR(:,:,it(j));
                        Dragforcethet(:,:,j)=Dragforcethet(:,:,j)+-m_e*dragfreqsrc.*obj.fluidUTHET(:,:,it(j));
                        Dragforcez(:,:,j)=Dragforcez(:,:,j)+-m_e*dragfreqsrc.*obj.fluidUZ(:,:,it(j));
                    end
                    
                    % Time derivative
                    cdt=(obj.t2d(it(j)+1)-obj.t2d(it(j)-1));
                    durdt(:,:,j)=m_e*(obj.fluidUR(:,:,it(j)+1)-obj.fluidUR(:,:,it(j)-1))/cdt;
                    duthetdt(:,:,j)=m_e*(obj.fluidUTHET(:,:,it(j)+1)-obj.fluidUTHET(:,:,it(j)-1))/cdt;
                    duzdt(:,:,j)=m_e*(obj.fluidUZ(:,:,it(j)+1)-obj.fluidUZ(:,:,it(j)-1))/cdt;
                end
                if(~fdens)
                    Forces.Eforcer=Eforcer;
                    Forces.Eforcez=Eforcez;
                    Forces.Bforcer=Bforcer;
                    Forces.Bforcethet=Bforcethet;
                    Forces.Bforcez=Bforcez;
                    Forces.inertforcer=inertforcer;
                    Forces.inertforcethet=inertforcethet;
                    Forces.inertforcez=inertforcez;
                    Forces.Pressforcer=Pressforcer;
                    Forces.Pressforcethet=Pressforcethet;
                    Forces.Pressforcez=Pressforcez;
                    Forces.durdt=durdt;
                    Forces.duthetdt=duthetdt;
                    Forces.duzdt=duzdt;
                    Forces.Dragforcer=Dragforcer;
                    Forces.Dragforcethet=Dragforcethet;
                    Forces.Dragforcez=Dragforcez;
                else
                    % multiply by density to have force density
                    Forces.Eforcer=Eforcer.*n;
                    Forces.Eforcez=Eforcez.*n;
                    Forces.Bforcer=Bforcer.*n;
                    Forces.Bforcethet=Bforcethet.*n;
                    Forces.Bforcez=Bforcez.*n;
                    Forces.inertforcer=inertforcer.*n;
                    Forces.inertforcethet=inertforcethet.*n;
                    Forces.inertforcez=inertforcez.*n;
                    Forces.Pressforcer=Pressforcer.*n;
                    Forces.Pressforcethet=Pressforcethet.*n;
                    Forces.Pressforcez=Pressforcez.*n;
                    Forces.durdt=durdt.*n;
                    Forces.duthetdt=duthetdt.*n;
                    Forces.duzdt=duzdt.*n;
                    Forces.Dragforcer=Dragforcer.*n;
                    Forces.Dragforcethet=Dragforcethet.*n;
                    Forces.Dragforcez=Dragforcez.*n;
                    
                end
            end
        end
        
        function [lr,rb,lz,zb]= clouddims(obj,it,zpos,fracn)
            % clouddims return the cloud axial and radial limit at time it
            % and axial position zpos
            % fracn defines the fraction of the maximum density below which
            % we consider to have a vacuum
            if nargin<4
                fracn=0.1;
            end
            % get the density
            n=obj.N(:,:,it);
            lr=cell(1,length(it));
            lz=lr;
            rb=lr;
            zb=rb;
            for i=1:size(n,3)
                nthresh=fracn*max(max(n(:,:,i)));
                % find the points outside of the cloud
                outside=find(n(:,zpos,i)<nthresh);
                gap=diff(outside);
                k=1;
                for j=1:length(gap)
                    if(gap(j)>2)
                        rmpos=outside(j);
                        rppos=outside(j+1);
                        lr{i}(k)=obj.rgrid(rppos-1)-obj.rgrid(rmpos+1);
                        rb{i}(:,k)=[max(rmpos+1,1) min(rppos-1,sum(obj.nnr))];
                        k=k+1;
                    end
                end
                maxgap=2;
                k=1;
                for I=rmpos+1:rppos-1
                    outside=find(n(I,:,i)<nthresh);
                    zgap=diff(outside);
                    for j=1:length(zgap)
                        if(zgap(j)>maxgap)
                            maxgap=zgap(j);
                            zmpos=outside(j);
                            zppos=outside(j+1);
                            lz{i}(k)=obj.zgrid(zppos-1)-obj.zgrid(zmpos+1);
                            zb{i}(:,k)=[max(zmpos+1,1) min(zppos-1,obj.nz)];
                            k=k+1;
                        end
                    end
                end
                
            end
        end
        
        %------------------------------------------
        %  Functions for plotting evolving quantities
        
        function displaysplbound(obj,ax,rescale)
            %displaysplbound display on axis ax the boundary of the
            %simulation domain and the Dirichlet and Neumann walls defined
            %with spline curves
            if nargin<2
                ax=gca;
            end
            if nargin<3
                rescale=1;
            end
            hold on
            for i=1:obj.spl_bound.nbsplines
                knots=obj.spl_bound.boundary(i).knots(1:end);
                coeffs=obj.spl_bound.boundary(i).coefs'*rescale;
                pp=spmak(knots,coeffs);
                sizec=size(coeffs,2);
                order=length(knots)-sizec;
                s=linspace(knots(order),knots(sizec+1),1000);
                fittedpos=fnval(pp,s);
                plot(fittedpos(1,:),fittedpos(2,:),'-')
                plot(coeffs(1,:),coeffs(2,:),'rx','markersize',14)
            end
        end
        
        function displayraddim(obj,it,zpos,fracn)
            %displayraddim display the evolution of the radial dimension of the cloud in
            %time to find if the cloud size get below a critical radial
            %size at which the ionisation is not sufficient to compensate
            %the losses
            % also plot the well radial dimensions in time
            if nargin<3
                zpos=floor(length(obj.zgrid)/2);
            end
            if nargin<4
                fracn=0.1;
            end
            [lr,rb,lz,zb]=obj.clouddims(it,zpos,fracn);
            
            t=obj.t2d(it);
            Lr=zeros(size(lr));
            er=obj.Er(:,:,it);
            r_min=Lr;
            r_minpred=r_min;
            well_r=Lr;
            nb=Lr;
            for i=1:length(lr)
                if ~isempty(lr{i}) && ~isempty(lz{i})
                    [Lr(i),id]=max(lr{i});
                    rm=rb{i}(1,id);
                    rp=rb{i}(2,id);
                    
                    nb(i)=mean(obj.N(rm:rp,zpos,it(i)));
                    
                    Lp=min(lz{i});
                    Lm=mean(lz{i});
                    
                    rpos=rm:rp;
                    vperp=-er(rpos,zpos,i)./obj.Bz(zpos,rpos)';
                    Ek=0.5*obj.me*vperp.^2/obj.qe;
                    sigio=obj.sigio(Ek);
                    sigd=obj.sigmela(Ek)+obj.sigmio(Ek)+sigio;
                    omegap2=obj.qe^2*obj.N(rpos,zpos,it(i))/obj.eps_0/obj.me;
                    omegac2=(obj.qe*obj.Bz(zpos,rpos)'/obj.me).^2;
                    ur=er(rpos,zpos,i)*obj.qe./((omegap2-omegac2)*obj.me).*sigd.*vperp*obj.neutcol.neutdens;
                    r_minpred(i)=mean(obj.N(rp,zpos,it(i))*Lp*ur./(nb(i)*obj.neutcol.neutdens*sigio.*vperp*Lm));%mean(1./(-1/obj.rgrid(rm)+obj.neutcol.neutdens*sigio.*vperp./ur*(Lm/Lp)*nb(i)/obj.N(rp,zpos,it(i))));
                    rpos=rp;
                    vperp=-er(rpos,zpos,i)./obj.Bz(zpos,rpos)';
                    Ek=0.5*obj.me*vperp.^2/obj.qe;
                    sigio=obj.sigio(Ek);
                    ur=obj.fluidUR(rpos,zpos,it(i));
                    r_min(i)=max(obj.N(rp,zpos,it(i))*Lp*ur/(nb(i)*obj.neutcol.neutdens*sigio.*vperp*Lm),0);%max(mean(1./(-1/obj.rgrid(rm)+obj.neutcol.neutdens*sigio.*vperp./ur*(Lm/Lp)*nb(i)/obj.N(rp,zpos,it(i)))),0);
                    
                    nb(i)=nb(i)*Lm*2*pi*obj.rgrid(rm)*Lr(i);
                else
                    Lr(i)=NaN;
                    r_min(i)=NaN;
                    r_minpred(i)=NaN;
                end
                potwell=obj.PotentialWell(it(i))';
                outside=find(isnan(potwell(:,zpos)));
                gap=diff(outside);
                for j=1:length(gap)
                    if(gap(j)>2)
                        rmpos=outside(j)+1;
                        rppos=outside(j+1)-1;
                        well_r(i)=obj.rgrid(rppos)-obj.rgrid(rmpos);
                    end
                end
                
                
                
            end
            
            f=figure('Name', sprintf('%s rlims B=%f phi=%f',obj.name,obj.B0, obj.phinorm*(obj.potout-obj.potinn)));
            plot(t,Lr,'displayname','\Deltar_{cloud}','linewidth',1.3)
            hold on
            plot(t,r_min,'displayname','\Deltar_{min} (u_r simu)','linewidth',1.3)
            plot(t,r_minpred,'displayname','\Deltar_{min} (u_r pred)','linewidth',1.3)
            plot(t,well_r,'displayname','\Deltar_{well}','linewidth',1.3)
            ylabel('\Delta r [m]')
            yyaxis right
            plot(t,nb,'--','displayname','N')
            legend('location','eastoutside')
            xlabel('t [s]')
            ylabel('N')
            set(gca,'fontsize',12)
            
            yyaxis left
            ylimits=ylim;
            %ylim([ylimits(1) 1.1*max(Lr)])
            title(sprintf('cloud radial limits at z=%1.2e[m]',obj.zgrid(zpos)))
            obj.savegraph(f,sprintf('%s/%s_%d_rlims',obj.folder,obj.name,zpos),[15 10]);
        end
        
        function displaypsi(obj,deltat)
            %% plot the initial and final radial profile at position z=0 and show the normalized enveloppe function Psi
            % relevant for Davidson annular distribution function
            f=figure('Name', sprintf('%s Psi',obj.name));
            f.Name= sprintf('%s Psi',obj.name);
            zpos=floor(length(obj.zgrid)/2);
            tinit=1;
            tend=length(obj.t2d);
            if iscell(deltat)
                deltat=cell2mat(deltat);
            end
            if(obj.R.nt<2)
                h0=obj.H0;
                p0=obj.P0;
            else
                h0=mean(H(obj,{1:obj.VR.nparts,obj.VR.nt,false}));
                p0=mean(P(obj,{1:obj.VR.nparts,obj.VR.nt,false}));
            end
            lw=1.5;
            Mirrorratio=(obj.Rcurv-1)/(obj.Rcurv+1);
            locpot=mean(obj.pot(:,zpos,tend-deltat:tend),3);
            psi=1+obj.qe*locpot(:)/h0-1/(2*obj.me*h0)*(p0./obj.rgrid+obj.qe*0.5*obj.B0.*(obj.rgrid-obj.width/pi*Mirrorratio*cos(2*pi*obj.zgrid(zpos)/obj.width)*besseli(1,2*pi*obj.rgrid/obj.width))).^2;
            locdens=mean(obj.N(:,zpos,tend-deltat:tend),3);
            [maxn,In]=max(locdens);%M.N(:,zpos,tinit));
            plot(obj.rgrid,obj.N(:,zpos,tinit),'bx-','DisplayName',sprintf('t=%1.2f[ns]',obj.t2d(tinit)*1e9),'linewidth',lw)
            hold on
            plot(obj.rgrid,locdens,'rx-','DisplayName',sprintf('t=[%1.2f-%1.2f] [ns] averaged',obj.t2d(tend-deltat)*1e9,obj.t2d(tend)*1e9),'linewidth',lw)
            plot(obj.rgrid(In-2:end),1./obj.rgrid(In-2:end)*maxn*obj.rgrid(In),'DisplayName','N=c*1/r','linewidth',lw)
            plot(obj.rgrid(In-2:end),1./obj.rgrid(In-2:end).^2*maxn*obj.rgrid(In)^2,'DisplayName','N=c*1/r^2','linewidth',lw)
            plot(obj.rgrid(In-2:end),1./obj.rgrid(In-2:end).^4*maxn*obj.rgrid(In)^4,'DisplayName','N=c*1/r^4','linewidth',lw)
            xlabel('r [m]')
            ylabel('n_e [m^{-3}]')
            I=find(psi>0);
            if (length(I)>1)
                I=[I(1)-2; I(1)-1; I; I(end)+1; I(end)+2];
            else
                I=obj.nnr(1):length(psi);
            end
            rq=linspace(obj.rgrid(max(I(1),1)),obj.rgrid(I(end)),500);
            psiinterp=interp1(obj.rgrid(I),psi(I),rq,'pchip');
            zeroindices=find(diff(psiinterp>=0),2);
            maxpsiinterp=max(psiinterp);
            plot(rq,maxn*psiinterp/abs(maxpsiinterp),'Displayname','normalized \Psi [a.u.]','linewidth',lw)
            ylim([0 inf])
            for i=1:length(zeroindices)
                border=plot([rq(zeroindices(i)) rq(zeroindices(i))],[0 obj.rgrid(In)/obj.rgrid(In-2)*maxn],'k--','linewidth',lw);
                set(get(get(border,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            legend
            xlim([0 0.02])
            grid on
            title(sprintf('Radial density profile at z=%1.2e[m]',obj.zgrid(zpos)))
            obj.savegraph(f,sprintf('%sPsi',obj.name),[15 10]);
        end
        
        function f=displayrprofile(obj,t,zpos,init)
            %% plot the initial and final radial profile at the axial center of the simulation space
            % also plot the azimuthal fluid rotation frequency profile
            % t: time index considered
            % zpos: axial position index
            % init: initial time considered for comparison
            f=figure('Name', sprintf('%s Prof',obj.name));
            if nargin < 3 || length(zpos)<1
                zpos=floor(length(obj.zgrid)/2);
            end
            if nargin<4
                init=false;
            end
            if(iscell(t))
                t=cell2mat(t);
            end
            lw=1.5;
            if init
                tinit=t(1);
                t=t(2:end);
            end
            locdens=mean(obj.N(:,zpos,t),3);
            
            %inverse of radius
            Rinv=1./obj.rgrid;
            Rinv(isnan(Rinv))=0;
            %azimuthal velocity and azimuthal rotation frequency in m/s and
            %1/s
            vthet=mean(obj.fluidUTHET(:,zpos,t),3);
            omegare=(vthet.*Rinv);
            % plot the initial density
            if(init)
                plot(obj.rgrid,obj.N(:,zpos,tinit),'bx-','DisplayName',sprintf('t=%1.2f[ns]',obj.t2d(tinit)*1e9),'linewidth',lw)
            end
            hold on
            %plot the time averaged current density
            plot(obj.rgrid,locdens,'rx-','DisplayName',sprintf('t=[%1.2f-%1.2f] [ns] averaged',obj.t2d(t(1))*1e9,obj.t2d(t(end))*1e9),'linewidth',lw)
            xlabel('r [m]')
            ylabel('n_e [m^{-3}]')
            legend('location','Northwest')
            % limit the axis to the simulation domain
            if obj.conformgeom
                xlim([obj.rgrid(1) obj.rgrid(sum(obj.nnr(1:2)))])
            else
                xlim([obj.rgrid(1) obj.rgrid(end)])
            end
            grid on
            ylimits=ylim();
            % plot the metallic walls for a constant radius coaxial
            % configuration
            if obj.conformgeom
                plot(obj.rgrid(1)*[1 1],ylimits,'k--')
                plot(obj.rgrid(end)*[1 1],ylimits,'k--')
            else
                plot(obj.r_a*[1 1],ylimits,'k--')
                if obj.walltype==0
                    plot(obj.r_b*[1 1],ylimits,'k--')
                elseif obj.walltype==1
                    rmax=obj.r_0-obj.r_r*sqrt(1-(obj.zgrid(zpos)-obj.z_0)^2/obj.z_r^2);
                    plot(rmax*[1 1],ylimits,'k--')
                end
            end
            
            yyaxis right
            % plot the azimuthal fluid rotation frequency profile
            plot(obj.rgrid,omegare,'DisplayName',sprintf('<\\omega_{re}> t=[%1.2f-%1.2f] [ns] averaged',obj.t2d(t(1))*1e9,obj.t2d(t(end))*1e9),'linewidth',lw)
            ylabel('\omega_{re} [1/s]')
            
            title(sprintf('Radial density profile at z=%1.2e[m]',obj.zgrid(zpos)))
            obj.savegraph(f,sprintf('%srProf',obj.name),[15 10]);
        end
        
        function displayenergy(obj)
            %% Plot the time evolution of the system energy and number of simulated macro particles
            tmin=1;
            tmax=length(obj.ekin);
            f=figure('Name', sprintf('%s Energy',obj.name));
            subplot(2,1,1)
            plot(obj.t0d(tmin:tmax),obj.ekin(tmin:tmax),'o-',...
                obj.t0d(tmin:tmax),obj.epot(tmin:tmax),'d-',...
                obj.t0d(tmin:tmax),obj.etot(tmin:tmax),'h-',...
                obj.t0d(tmin:tmax),obj.etot0(tmin:tmax),'h-',...
                obj.t0d(tmin:tmax),obj.eerr(tmin:tmax),'x--')
            %obj.t0d(tmin:tmax),obj.ekin(tmin:tmax)-obj.epot(tmin:tmax),'--',
            legend('E_{kin}', 'E_{pot}', 'E_{tot}','E_{ref}','E_{err}')
            xlabel('Time [s]')
            ylabel('Energies [J]')
            grid on
            xlimits=xlim();
            
            subplot(2,1,2)
            try
                semilogy(obj.t0d(tmin:tmax),abs(obj.eerr(tmin:tmax)./obj.etot0(tmin:tmax)),'h-')
            catch
                semilogy(obj.t0d(tmin:tmax),abs(obj.eerr(tmin:tmax)/obj.etot(2)),'h-')
            end
            hold on
            xlabel('t [s]')
            ylabel('E_{err}/E_{tot}')
            xlim(xlimits)
            grid on
            try
                yyaxis right
                plot(obj.t0d(tmin:tmax),abs(obj.npart(tmin:tmax)./obj.npart(1)*100),'d--')
                ylabel('Nparts %')
                %ylim([0 110])
            catch
            end
            ylimits=ylim;
            for i=1:length(obj.restarttimes)
                plot(obj.restarttimes(i)*[1 1],ylimits,'k--')
            end
            
            obj.savegraph(f,sprintf('%s/%sEnergy',obj.folder,obj.name));
        end
        
        function f=displaycharge(obj,f,linelegend)
            %% Plot the time evolution of the system charge of electrons
            % f: figure handle if you want to stack several such curves
            % on the same figure
            % linelegend: legend for this charge evolution
            tmin=1;
            tmax=length(obj.ekin);
            if nargin<2
                f=figure('Name', sprintf('%s Charge',obj.name));
            end
            if nargin<3
                linelegend='';
            end
            ax=f.CurrentAxes;
            if isempty(ax)
                ax=axes(f);
            end
            try
                plot(ax,obj.t0d(tmin:tmax),abs(obj.npart(tmin:tmax)*obj.qsim),'linewidth',1.5,'displayname',linelegend)
                hold on
                ylabel(ax,'Total charge [C]')
                xlabel(ax,'t [s]')
                grid on
                if(nargin>2)
                    legend
                end
                set(ax,'fontsize',12)
            catch
            end
            if nargin < 2
                obj.savegraph(f,sprintf('%s/%scharge',obj.folder,obj.name));
            end
        end
        
        function displaySimParticles(obj)
            %% Plot the time evolution of the number of simulated markers in the main specie
            f=figure('Name', sprintf('%s Trapped particles',obj.name));
            plot(obj.t0d,obj.npart,'linewidth',1.5)
            xlabel('t [s]')
            ylabel('N particles')
            set(gca,'fontsize',12)
            obj.savegraph(f,sprintf('%s/%sntrapped',obj.folder,obj.name),[10 12]);
        end
        
        function displayLarmorRad(obj,time2d)
            if nargin<2
                time2d=length(obj.t2d);
            end
            % Plot the larmor radius for created particles with low energy
            % the larmor radius is calculated by considering that the
            % initial perpendicular velocity \approx the ExB velocity
            if time2d>0
                Er=obj.Er(:,:,time2d);
                Ez=obj.Ez(:,:,time2d);
            else
                Er=obj.Erxt(:,:,1);
                Ez=obj.Ezxt(:,:,1);
            end
            
            rl=abs(obj.me/obj.qe*(-Er.*obj.Bz'+Ez.*obj.Br')./(obj.B.^3)');
            figure
            rl(obj.geomweight(:,:,1)<0)=0;
            contourf(obj.zgrid,obj.rgrid,rl)
            hold on
            contour(obj.zgrid,obj.rgrid,obj.geomweight(:,:,1),[0 0],'r-','linewidth',3)
            if time2d>0
                n=obj.N(:,:,time2d);
                maxN=max(n (:));
                n=n/maxN*mean(rl(:));
                contour(obj.zgrid,obj.rgrid,n,linspace(0,1,6)*mean(rl(:)),'r:','linewidth',3)
            end
            c=colorbar;
            xlabel('z [m]')
            ylabel('r [m]')
            c.Label.String='r_L [m]';
        end
        
        function displayHP(obj,tstart)
            % Plot the histogramm of the total energy and canonical angular momentum at time tstart and
            % end time of the simulation over the full simulation space for the main specie
            if(iscell(tstart))
                tstart=cell2mat(tstart);
            end
            if(obj.R.nt>=2)
                tstart=obj.R.nt;
                f=figure('Name', sprintf('%s HP',obj.name));
                legtext=sprintf("t=%2.1f - %2.1f [ns]",obj.tpart(tstart)*1e9,obj.tpart(end)*1e9);
                subplot(1,2,1)
                partsmax=min(obj.nbparts(end),obj.R.nparts);
                Hloc=H(obj,{1:obj.nbparts(1),1,false});
                h1=histogram(Hloc,20,'BinLimits',[min(Hloc(:)) max(Hloc(:))],'DisplayName',sprintf("t=%2.3d [ns]",obj.tpart(1)*1e9));
                hold on
                Hloc=H(obj,{1:partsmax,obj.R.nt,false});
                %,'Binwidth',h1.BinWidth
                h1=histogram(Hloc,20,'BinLimits',[min(Hloc(:)) max(Hloc(:))],'DisplayName',legtext);
                ylabel('counts')
                xlabel('H [J]')
                legend
                
                subplot(1,2,2)
                Ploc=P(obj,{1:obj.nbparts(1),1,false});
                h2=histogram(Ploc,50,'BinLimits',[min(Ploc(:)) max(Ploc(:))],'DisplayName',sprintf("t=%2.3d [ns]",obj.tpart(1)*1e9));
                hold on
                Ploc=P(obj,{1:partsmax,obj.R.nt,false});
                histogram(Ploc,50,'BinLimits',[min(Ploc(:)) max(Ploc(:))],'DisplayName',legtext);
                ylabel('counts')
                xlabel('P [kg\cdotm^2\cdots^{-1}]')
                %clear P
                %clear H
                legend
                %xlim([0.95*h2.BinLimits(1) 1.05*h2.BinLimits(2)])
                obj.savegraph(f,sprintf('%s/%sParts_HP',obj.folder,obj.name));
            end
        end
        
        function displayaveragetemp(obj)
            % Computes and show the particles average temperature as a function of time
            f=figure('Name',sprintf('%s potinn=%f part temperature',obj.name,obj.potinn*obj.phinorm));
            vr2=obj.VR(:,:,false);
            vr2=mean(vr2.^2,1)-mean(vr2,1).^2;
            vz2=obj.VZ(:,:,false);
            vz2=mean(vz2.^2,1)-mean(vz2,1).^2;
            vthet2=obj.VTHET(:,:,false);
            vthet2=mean(vthet2.^2,1)-mean(vthet2,1).^2;
            plot(obj.tpart,0.5*obj.me*vr2/obj.qe,'displayname','T_r')
            hold on
            plot(obj.tpart,0.5*obj.me*vz2/obj.qe,'displayname','T_z')
            plot(obj.tpart,0.5*obj.me*vthet2/obj.qe,'displayname','T_{thet}')
            xlabel('time [s]')
            ylabel('T [eV]')
            title(sprintf('\\phi_a=%.1f kV \\phi_b=%.1f kV R=%.1f',obj.potinn*obj.phinorm/1e3,obj.potout*obj.phinorm/1e3,obj.Rcurv))
            legend
            grid
            obj.savegraph(f,sprintf('%s/%s_partstemp',obj.folder,obj.name));
            
        end
        
        function displayCurrentsevol(obj,timesteps)
            % Computes and display the time evolution of the outgoing currents on each domain boundary
            % at timesteps timesteps
            if nargin<2
                timesteps=1:length(obj.t2d);
            end
            currents=obj.OutCurrents(timesteps);
            f=figure('Name',sprintf('%s Currents',obj.name));
            if(obj.B(1,1)>obj.B(end,1))
                lname='HFS';
                rname='LFS';
            else
                lname='LFS';
                rname='HFS';
            end
            plot(obj.t2d(timesteps),currents(1,:),'Displayname',lname,'linewidth',1.8);
            hold on
            plot(obj.t2d(timesteps),currents(2,:),'Displayname',rname,'linewidth',1.8);
            plot(obj.t2d(timesteps),currents(3,:),'Displayname','outer cylinder','linewidth',1.8);
            plot(obj.t2d(timesteps),currents(4,:),'Displayname','inner cylinder','linewidth',1.8);
            if size(currents,1)>=5
                plot(obj.t2d(timesteps),currents(5,:),'Displayname','ellipse','linewidth',1.8);
            end
            legend('location','Northeast')
            xlabel('time [s]')
            ylabel('I [A]')
            grid on
            set(gca,'fontsize',12)
            title(sprintf('\\phi_b-\\phi_a=%.2g kV, R=%.1f',(obj.potout-obj.potinn)*obj.phinorm/1e3,obj.Rcurv))
            obj.savegraph(f,sprintf('%s/%s_outCurrents',obj.folder,obj.name),[16 12]);
        end
        
        function displayChargeLossevol(obj,timesteps,toptitle,scalet,dens)
            % Computes and display the time evolution of the outgoing currents on each domain boundary
            % at time obj.t2d(timesteps)
            %scalet=true scales the time by the ellastic collision
            %frequency
            %dens = true plot the time evolution of the maximum electron
            %density in the simulation domain otherwise plot the total
            %number of electrons in the domain
            if nargin<2
                timesteps=1:length(obj.t2d);
            end
            if nargin<4
                scalet=true;
            end
            if nargin <5
                dens=true;
            end
            
            if scalet
                if obj.neutcol.present
                    vexb0=(obj.Ez(:,:,1).*obj.Br'-obj.Er(:,:,1).*obj.Bz')./(obj.B'.^2);
                    vexb0(obj.geomweight(:,:,1)<=0)=0;
                    E=0.5*obj.msim/obj.weight*mean(abs(vexb0(:)))^2/obj.qe;
                    taucol=1/(obj.neutcol.neutdens*mean(abs(vexb0(:)))*(obj.sigio(E)+obj.sigmela(E)+obj.sigmio(E)));
                    try
                        Sio_S=1e17*(obj.neutcol.neutdens*mean(abs(vexb0(:)))*obj.sigio(E))/(obj.maxwellsrce.frequency*obj.weight/(pi*(obj.maxwellsrce.rlim(2)^2-obj.maxwellsrce.rlim(1)^2)*diff(obj.maxwellsrce.zlim)))
                    catch
                    end
                    tlabel='t/\tau_d [-]';
                else
                    taucol=2*pi/obj.omece;
                    tlabel='t/\tau_ce [-]';
                end
            else
                taucol=1e-9;
                tlabel='t [ns]';
            end
            
            if dens
                N=obj.N(:,:,timesteps);
                geomw=obj.geomweight(:,:,1);
                geomw(geomw<0)=0;
                geomw(geomw>0)=1;
                N=N.*geomw;
                
                nmax=squeeze(max(max(N,[],1),[],2));
                tn=(obj.t2d(timesteps));
                nlabel='n_{e,max} [m^{-3}]';
                ndlabel='n_{e,max}';
            else
                t0dst=find(obj.t0d>=obj.t2d(timesteps(1)),1,'first');
                t0dlst=find(obj.t0d<=obj.t2d(timesteps(end)),1,'last');
                tn=obj.t0d(t0dst:t0dlst);
                nmax=obj.npart(t0dst:t0dlst)*obj.weight;
                nlabel='Nb e^-';
                ndlabel='Nb e^-';
            end
            
            
            [currents,pos]=obj.OutCurrents(timesteps);
            P=obj.neutcol.neutdens*obj.kb*300/100;% pressure at room temperature in mbar
            currents=currents/P;
            f=figure('Name',sprintf('%s Charges',obj.name));
            % Plot the evolution of nb of particles
            yyaxis right
            p=plot(tn/taucol,nmax,'b-.','linewidth',1.8,'Displayname',ndlabel);
            ylabel(nlabel)
            ax=gca;
            ax.YAxis(2).Color=p.Color;
            ylim([0 inf])
            
            if(obj.B(1,1)>obj.B(end,1))
                lname='HFS';
                rname='LFS';
            else
                lname='LFS';
                rname='HFS';
            end
            
            yyaxis left
            mincurr=max(currents(:))*5e-3;
            if (max(currents(1,:)>mincurr))
                plot(obj.t2d(timesteps)/taucol,currents(1,:),'r:','Displayname',lname,'linewidth',1.8);
            end
            hold on
            if (max(currents(2,:)>mincurr))
                plot(obj.t2d(timesteps)/taucol,currents(2,:),'r--','Displayname',rname,'linewidth',1.8);
            end
            if (max(currents(3,:)>mincurr))
                plot(obj.t2d(timesteps)/taucol,currents(3,:),'r-','Displayname','outer cylinder','linewidth',1.8);
            end
            if (max(currents(4,:)>mincurr))
                plot(obj.t2d(timesteps)/taucol,currents(4,:),'Displayname','inner cylinder','linewidth',1.8);
            end
            if (size(currents,1)>=5 && max(currents(5,:)>mincurr))
                plot(obj.t2d(timesteps)/taucol,currents(5,:),'r-','Displayname','ellipse','linewidth',1.8);
            end
            xlabel(tlabel)
            ylabel('I/p_n [A/mbar]')
            grid on
            set(gca,'fontsize',12)
            ax.YAxis(1).Color='red';
            
            legend('Orientation','horizontal','location','south','numcolumns',3)
            
            if nargin <3
                title(sprintf('\\phi_b-\\phi_a=%.2g kV, B=%f T',(obj.potout-obj.potinn)*obj.phinorm/1e3,max(obj.B(:))))
            elseif ~isempty(toptitle)
                title(toptitle)
            end
            obj.savegraph(f,sprintf('%s/%s_ChargeEvol%i%i',obj.folder,obj.name,scalet,dens),[16 14]);
        end
        
        
        function displaytotcurrevol_geom(obj,timesteps,toptitle,scalet,dens,subdiv,nmean)
            % Computes and display the time evolution of the outgoing
            % currents at time obj.t2d(timesteps)
            %scalet=true scales the time by the ellastic collision
            %frequency
            %dens = true plot the time evolution of the maximum electron
            %density in the simulation domain otherwise plot the total
            %number of electrons in the domain
            % also plot in a subplot the color coded boundary corresponding
            % to each current
            if nargin<2
                timesteps=1:length(obj.t2d);
            end
            if nargin<3
                scalet=true;
            end
            if nargin <4
                dens=true;
            end
            if nargin<5
                subdiv=1;
            end
            if nargin<6
                nmean=1;
            end
            
            if scalet
                if obj.neutcol.present
                    vexb0=(obj.Ezxt(:,:,1).*obj.Br'-obj.Erxt(:,:,1).*obj.Bz')./(obj.B'.^2);
                    vexb0(obj.geomweight(:,:,1)<=0)=0;
                    potwell=obj.PotentialWell(0);
                    vexb0(isnan(potwell))=NaN;
                    vexb0=mean(abs(vexb0(:)),'omitnan');
                    E=0.5*obj.msim/obj.weight*vexb0^2/obj.qe;
                    taucol=1/(obj.neutcol.neutdens*vexb0*(obj.sigio(E)+obj.sigmela(E)+obj.sigmio(E)));
                    try
                        Sio_S=1e17*(obj.neutcol.neutdens*vexb0*obj.sigio(E))/(obj.maxwellsrce.frequency*obj.weight/(pi*(obj.maxwellsrce.rlim(2)^2-obj.maxwellsrce.rlim(1)^2)*diff(obj.maxwellsrce.zlim)))
                    catch
                    end
                    tlabel='t/\tau_d [-]';
                else
                    taucol=2*pi/obj.omece;
                    tlabel='t/\tau_ce [-]';
                end
            else
                taucol=1e-9;
                tlabel='t [ns]';
            end
            
            if dens
                N=obj.N(:,:,timesteps);
                geomw=obj.geomweight(:,:,1);
                geomw(geomw<0)=0;
                geomw(geomw>0)=1;
                N=N.*geomw;
                %[~,idl]=max(N(:,:,end),[],'all','linear');
                %[ir,iz]=ind2sub(size(geomw),idl);
                %nmax=squeeze(max(max(N,[],1),[],2));
                tn=(obj.t2d(timesteps));
                nrhalf=find(obj.rgrid>0.07,1,'first');
                if(nrhalf<length(obj.rgrid))
                    nmax=zeros(2,length(tn));
                    
                    nmax(1,:)=squeeze(max(max(N(1:nrhalf,:,:),[],1),[],2));
                    nmax(2,:)=squeeze(max(max(N(nrhalf+1:end,:,:),[],1),[],2));
                else
                    nmax(:)=squeeze(max(max(N(:,:,:),[],1),[],2));
                end
                
                
                nlabel='n_{e,max} [m^{-3}]';
                ndlabel='n_{e,max}';
            else
                t0dst=find(obj.t0d>=obj.t2d(timesteps(1)),1,'first');
                t0dlst=find(obj.t0d<=obj.t2d(timesteps(end)),1,'last');
                tn=obj.t0d(t0dst:t0dlst);
                nmax=obj.npart(t0dst:t0dlst)*obj.weight;
                nlabel='Nb e^-';
                ndlabel='Nb e^-';
            end
            
            
            [currents,pos]=obj.OutCurrents(timesteps,subdiv);
            P=obj.neutcol.neutdens*obj.kb*300/100;% pressure at room temperature in mbar
            currents=currents/P;
            f=figure('Name',sprintf('%s Charges',obj.name));
            tiledlayout(2,1)
            nexttile
            % Plot the evolution of nb of particles
            yyaxis right
            if size(nmax,2)>1
                nmax=nmax';
            end
            for i=1:size(nmax,2)
                p=plot(tn/taucol,nmax(:,i),'b','linewidth',1.8,'Displayname',sprintf('%s, %d',ndlabel,i));
                hold on
            end
            ylabel(nlabel)
            axl=gca;
            axl.YAxis(2).Color=p.Color;
            ylim([0 inf])
            
            if(obj.B(1,1)>obj.B(end,1))
                lname='HFS';
                rname='LFS';
            else
                lname='LFS';
                rname='HFS';
            end
            
            yyaxis( 'left');
            map=colormap(lines);
            set(axl,'linestyleorder',{'-',':','--','*','+'},...
                'ColorOrder',map(2:7,:), 'NextPlot','replacechildren')
            p(1)=plot(axl,obj.t2d(timesteps)/taucol,movmean(currents(1,:),nmean),'Displayname',lname,'linewidth',1.8);
            hold on
            p(2)=plot(axl,obj.t2d(timesteps)/taucol,movmean(currents(2,:),nmean),'Displayname',rname,'linewidth',1.8);
            % Plot the currents
            for i=3:size(currents,1)
                p(i)=plot(axl,obj.t2d(timesteps)/taucol,movmean(currents(i,:),nmean),'Displayname',sprintf('border %i',i-2),'linewidth',1.8);
            end
            plot(axl,obj.t2d(timesteps)/taucol,movmean(sum(currents(:,:),1,'omitnan'),nmean),'k-','Displayname','total','linewidth',1.8);
            xlabel(tlabel)
            ylabel('I/p_n [A/mbar]')
            grid on
            set(gca,'fontsize',12)
            ax.YAxis(1).Color='black';
            
            %legend('Orientation','horizontal','location','north','numcolumns',3)
            
            if ~isempty(toptitle)
                title(toptitle)
            end
            
            ax2=nexttile;
            geomw=obj.geomweight(:,:,1);
            geomw(geomw<0)=0;
            geomw(geomw>0)=NaN;
            [c1,hContour]=contourf(ax2,obj.zgrid*1000,obj.rgrid*1000,geomw, [0 0]);
            hold on
            drawnow;
            grid on;
            
            for i=1:length(pos)
                plot(ax2,pos{i}(1,:)*1000,pos{i}(2,:)*1000,'linestyle',p(i+2).LineStyle,...
                    'color',p(i+2).Color,'marker',p(i+2).Marker,...
                    'displayname',sprintf('border %i',i),'linewidth',1.8)
                hold on
            end
            title('Domain')
            plot(ax2,ones(size(obj.rgrid))*obj.zgrid(1)*1000,obj.rgrid*1000,'linestyle',p(1).LineStyle,...
                'color',p(1).Color,'marker',p(1).Marker,...
                'displayname',lname,'linewidth',1.8)
            plot(ax2,ones(size(obj.rgrid))*obj.zgrid(end)*1000,obj.rgrid*1000,'linestyle',p(2).LineStyle,...
                'color',p(2).Color,'marker',p(2).Marker,...
                'displayname',rname,'linewidth',1.8)
            xlabel('z [mm]')
            ylabel('r [mm]')
            grid on
            set(gca,'fontsize',12)
            hFills=hContour.FacePrims;
            [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
            try
                hFills(1).ColorData = uint8([150;150;150;255]);
                for idx = 2 : numel(hFills)
                    hFills(idx).ColorData(4) = 0;   % default=255
                end
            catch
            end
            %legend('Orientation','horizontal','location','north','numcolumns',4)
            
            
            fprintf('mean total current: %f [A/mbar]\n',mean(sum(currents(:,max(1,size(currents,2)-30):end),1,'omitnan')));
            
            %if nargin <3
            %    sgtitle(sprintf('\\phi_b-\\phi_a=%.2g kV, B=%f T',(obj.potout-obj.potinn)*obj.phinorm/1e3,mean(obj.B(:))))
            %elseif ~isempty(toptitle)
            %    sgtitle(toptitle)
            %end
            if length(subdiv)>1
                obj.savegraph(f,sprintf('%s/%s_totIEvol%i%i_div%i',obj.folder,obj.name,scalet,dens,nmean),[16 14]);
            else
                obj.savegraph(f,sprintf('%s/%s_totIEvol%i%i_%i',obj.folder,obj.name,scalet,dens,nmean),[16 14]);
            end
            
        end
        
        function display1Dpotentialwell(obj,timestep,rpos)
            % Display the potential well along the magentic field line
            % passing by rgrid(rpos) at the center of the simulation space
            if iscell(timestep)
                timestep=cell2mat(timestep);
            end
            f=figure('Name',sprintf('%s 1D Potential well',obj.name));
            model=obj.potentialwellmodel(timestep);
            z=model.z;
            r=model.r;
            Pot=model.pot;
            rathet=model.rathet;
            if (mod(rpos, 1) ~= 0)
                [~,rpos]=min(abs(M.rgrid-rpos));
            end
            crpos=obj.rgrid(rpos);
            id=find(timestep==0);
            timestep(id)=[];
            n=obj.N(:,:,timestep);
            if(~isempty(timestep==0))
                N0=zeros(obj.N.nr+1,obj.N.nz+1);
                n=cat(3,n(:,:,1:id-1),N0,n(:,:,id:end));
            end
            n=mean(n,3);
            linepot=zeros(length(obj.zgrid),length(timestep));
            rathetpos=obj.rAthet(rpos,ceil(length(obj.zgrid)/2));
            F=scatteredInterpolant(z',rathet',Pot(:,1));
            for i=1:length(timestep)
                F=scatteredInterpolant(z',rathet',Pot(:,i));
                linepot(:,i)=F(obj.zgrid,rathetpos*ones(size(obj.zgrid)));
                %linepot(:,i)=griddata(z,rathet,pot(:,i),obj.zgrid,rathetpos);
            end
            linepot=mean(linepot,2);
            [Zinit,~]=meshgrid(obj.zgrid,obj.rAthet(:,1));
            n=griddata(Zinit,obj.rAthet,n,obj.zgrid,rathetpos);
            
            plot(obj.zgrid,linepot)
            ylabel('Potentiel [eV]')
            xlabel('z [m]')
            xlim([obj.zgrid(1) obj.zgrid(end)])
            hold(gca, 'on')
            yyaxis right
            plot(obj.zgrid,n)
            ylabel('n [m^{-3}]')
            if length(timestep)==1
                title(sprintf('Potential well t=%1.2f [ns] r=%1.2f [mm]',obj.t2d(timestep)*1e9,1e3*crpos))
            else
                title(sprintf('Potential well t=[%1.2f-%1.2f] [ns] r=%1.2f [mm]',obj.t2d(timestep(1))*1e9,obj.t2d(timestep(end))*1e9,1e3*crpos))
            end
            
            obj.savegraph(f,sprintf('%s/%s_well1Dr_%d',obj.folder,obj.name,rpos));
        end
        
        function displayVdistribRThetZ(obj,timestep, rpos, zpos)
            %displayVdistribRThetZ plot the velocity distribution function
            % in m/s
            %extracted from the markers at position window from rpos(1)
            %rpos(end) and zpos(1) to zpos(end)
            % and at time obj.tpart(timestep)
            %rpos and zpos are given as grid indices
            if(obj.R.nt>=2)
                if nargin<2
                    timesteppart=length(obj.tpart);
                else
                    timesteppart=timestep;
                end
                if nargin<3 || isempty(rpos)
                    rpos=1:length(obj.rgrid);
                    rspan=[obj.rgrid(1) obj.rgrid(end)];
                else
                    r=[obj.rgrid(1);(obj.rgrid(1:end-1)+obj.rgrid(2:end))*0.5;obj.rgrid(end)];
                    rspan=[r(rpos) r(rpos+1)];
                end
                if nargin<4 || isempty(zpos)
                    zpos=1:length(obj.zgrid);
                    zspan=[obj.zgrid(1) obj.zgrid(end)];
                else
                    z=[obj.zgrid(1);(obj.zgrid(1:end-1)+obj.zgrid(2:end))*0.5;obj.zgrid(end)];
                    zspan=[z(zpos) z(zpos+1)];
                end
                
                
                nbp=min(obj.nbparts(1),obj.R.nparts);
                R=obj.R(1:nbp,1,false);
                Z=obj.Z(1:nbp,1,false);
                Vr=obj.VR(1:nbp,1,false);
                Vz=obj.VZ(1:nbp,1,false);
                Vthet=obj.VTHET(1:nbp,1,false);
                ids=R>=rspan(1) & R<=rspan(2) & Z>=zspan(1) & Z<=zspan(2);
                
                Vr=Vr(ids);
                Vz=Vz(ids);
                Vthet=Vthet(ids);
                vTr=std(Vr,1);
                vTz=std(Vz,1);
                vTthet=std(Vthet,1);
                
                nbp=min(obj.nbparts(timesteppart),obj.R.nparts);
                Rend=obj.R(1:nbp,timesteppart,false);
                Zend=obj.Z(1:nbp,timesteppart,false);
                Vrend=obj.VR(1:nbp,timesteppart,false);
                Vzend=obj.VZ(1:nbp,timesteppart,false);
                Vthetend=obj.VTHET(1:nbp,timesteppart,false);
                ids=Rend>=rspan(1) & Rend<=rspan(2) & Zend>=zspan(1) & Zend<=zspan(2);
                nbtot=sum(ids)
                Vrend=Vrend(ids);
                Vzend=Vzend(ids);
                Vthetend=Vthetend(ids);
                vTrend=std(Vrend,1);
                vTzend=std(Vzend,1);
                vTthetend=std(Vthetend,1);
                
                binwidth=abs(max(Vrend)-min(Vrend))/sqrt(length(Vrend));
                f=figure('Name',sprintf("%s vrz distrib",obj.file));
                
                [~,time2did]=min(abs(obj.t2d-obj.tpart(timestep)));
                
                
                subplot(1,4,1);
                obj.dispV(Vr,Vrend,'V_r [m/s]',[1,timesteppart])
                
                [~,time2did]=min(abs(obj.t2d-obj.tpart(timestep)));
                if length(rpos)==1
                    vexb=-obj.Er(rpos,zpos,time2did)/obj.Bz(zpos,rpos)';
                    vexb=mean(vexb(:));
                    
                    if ~isempty(obj.neutcol.ela_cross_sec) % plot the radial drift velocity as nu_dE_r/(B\Omega_c)
                        vdr=obj.neutcol.neutdens*obj.sigmela(vexb^2*obj.me*0.5/obj.qe)*vexb*-obj.Er(rpos,zpos,time2did)...
                            ./(obj.B(zpos,rpos)'.*obj.B(zpos,rpos)'*obj.qe/obj.me);
                        vdr=mean(vdr(:));
                        ylimits=ylim;
                        plot(vdr*[1 1],ylimits,'k--','displayname',sprintf('V_{d,pred}=%1.2g [m/s]',vdr))
                    end
                end
                
                
                subplot(1,4,2);
                obj.dispV(Vthet,Vthetend,'V\theta [m/s]',[1,timesteppart])
                hold on
                drawnow
                ylimits=ylim;
                if length(rpos)==1
                    if ~isempty(obj.Erxt)
                        vexbext=-obj.Erxt(rpos,zpos)/obj.Bz(zpos,rpos)';
                        plot(vexbext*[1 1],ylimits,'k--','displayname',sprintf('V_{ExB,ext}=%1.2g [m/s]',vexbext))
                    end
                    
                    plot(vexb*[1 1],ylimits,'k-.','displayname',sprintf('V_{ExB,tot}=%1.2g [m/s]',vexb))
                end
                
                subplot(1,4,3);
                obj.dispV(Vz,Vzend,'Vz [m/s]',[1,timesteppart])
                
                subplot(1,4,4);
                obj.dispV(sqrt(Vr.^2+(Vthet).^2+Vz.^2),sqrt(Vrend.^2+(Vthetend).^2+Vzend.^2),'Vtot [m/s]',[1,timesteppart],'maxwell')
                
                
                sgtitle(sprintf('t=%1.2e[ns] r=[%2.1f, %2.1f] [mm] z=[%2.1f, %2.1f] [mm]',obj.tpart(timestep)*1e9, rspan*1e3, zspan*1e3))
                obj.savegraph(f,sprintf('%s/%sParts_V_RZ',obj.folder,obj.name),[25 14]);
                
                
            end
            
        end
        
        function displayEkin(obj,timestep, rpos, zpos)
            %displayEkin plot the kinetic energy distribution function in
            %eV
            %extracted from the markers at position window from rpos(1)
            %rpos(end) and zpos(1) to zpos(end)
            % and at time obj.tpart(timestep)
            %rpos and zpos are given as grid indices
            if(obj.R.nt>=2)
                if nargin<2
                    timesteppart=[1 length(obj.tpart)];
                else
                    if length(timestep)<2
                        timesteppart=[1 timestep];
                    else
                        timesteppart=[timestep(1) timestep(end)];
                    end
                end
                if nargin<3 || isempty(rpos)
                    rspan=[obj.rgrid(1) obj.rgrid(end)];
                else
                    r=[obj.rgrid(1);(obj.rgrid(1:end-1)+obj.rgrid(2:end))*0.5;obj.rgrid(end)];
                    rspan=[r(rpos) r(rpos+1)];
                end
                if nargin<4 || isempty(zpos)
                    zspan=[obj.zgrid(1) obj.zgrid(end)];
                else
                    z=[obj.zgrid(1);(obj.zgrid(1:end-1)+obj.zgrid(2:end))*0.5;obj.zgrid(end)];
                    zspan=[z(zpos) z(zpos+1)];
                end
                
                
                nbp=min(obj.nbparts(timesteppart(1)),obj.R.nparts);
                R=obj.R(1:nbp,timesteppart(1),false);
                Z=obj.Z(1:nbp,timesteppart(1),false);
                Ekin=obj.Ekin(1:nbp,timesteppart(1),false);
                ids=R>=rspan(1) & R<=rspan(2) & Z>=zspan(1) & Z<=zspan(2);
                Ekin=Ekin(ids)/obj.qe;
                
                nbp=min(obj.nbparts(timesteppart(2)),obj.R.nparts);
                Rend=obj.R(1:nbp,timesteppart(2),false);
                Zend=obj.Z(1:nbp,timesteppart(2),false);
                Ekinend=obj.Ekin(1:nbp,timesteppart(2),false);
                ids=Rend>=rspan(1) & Rend<=rspan(2) & Zend>=zspan(1) & Zend<=zspan(2);
                Ekinend=Ekinend(ids)/obj.qe;
                
                f=figure('Name',sprintf("%s E_k distrib",obj.file));
                
                obj.dispV(Ekin,Ekinend,'E_k',timesteppart,'maxwell')
                
                
                
                sgtitle(sprintf('dt=%1.2e[ns] r=[%2.1f, %2.1f] [mm] z=[%2.1f, %2.1f] [mm]',obj.dt*1e9, rspan*1e3, zspan*1e3))
                obj.savegraph(f,sprintf('%s/%sParts_E_kin',obj.folder,obj.name));
                
                
            end
            
        end
        
        function displayVdistribParPer(obj,timestep, rpos, zpos, gcs)
            %displayVdistribParPer plot the velocity distribution function
            % in m/s for the parallel and perpendicular velocity
            %extracted from the markers at position window from rpos(1)
            %rpos(end) and zpos(1) to zpos(end)
            % and at time obj.tpart(timestep)
            %rpos and zpos are given as grid indices
            % gcs define if you give the perpendicular velocity in the
            % guiding center frame or in the laboratory frame
            if(obj.R.nt>=2)
                if nargin<2
                    timesteppart=length(obj.tpart);
                else
                    timesteppart=timestep;
                end
                if nargin<3 || isempty(rpos)
                    rspan=[obj.rgrid(1) obj.rgrid(end)];
                else
                    r=[obj.rgrid(1);(obj.rgrid(1:end-1)+obj.rgrid(2:end))*0.5;obj.rgrid(end)];
                    rspan=[r(rpos) r(rpos+1)];
                end
                if nargin<4 || isempty(zpos)
                    zspan=[obj.zgrid(1) obj.zgrid(end)];
                else
                    z=[obj.zgrid(1);(obj.zgrid(1:end-1)+obj.zgrid(2:end))*0.5;obj.zgrid(end)];
                    zspan=[z(zpos) z(zpos+1)];
                end
                if nargin<5
                    gcs=false; % define if we look in the guiding center system
                end
                
                nbp=min(obj.nbparts(1),obj.R.nparts);
                R=obj.R(1:nbp,1,false);
                Z=obj.Z(1:nbp,1,false);
                ids=R>=rspan(1) & R<=rspan(2) & Z>=zspan(1) & Z<=zspan(2);
                Vperp=obj.Vperp(1:nbp,1,false,gcs);
                Vpar=obj.Vpar(1:nbp,1,false);
                Vperp=Vperp(ids);
                Vpar=Vpar(ids);
                
                nbp=min(obj.nbparts(timesteppart),obj.R.nparts);
                R=obj.R(1:nbp,timesteppart,false);
                Z=obj.Z(1:nbp,timesteppart,false);
                ids=R>=rspan(1) & R<=rspan(2) & Z>=zspan(1) & Z<=zspan(2);
                Vperpend=obj.Vperp(1:nbp,timesteppart,false,gcs);
                Vparend=obj.Vpar(1:nbp,timesteppart,false);
                Vperpend=Vperpend(ids);
                Vparend=Vparend(ids);
                
                %binwidth=abs(max(Vparend)-min(Vparend))/500;
                f=figure('Name',sprintf("%s v parper distrib",obj.file));
                subplot(1,2,1)
                if gcs
                    lgd='v_\perp gcs [m/s]';
                else
                    lgd='v_\perp [m/s]';
                end
                obj.dispV(Vperp,Vperpend,lgd,[1,timesteppart], 'maxwell')
                
                subplot(1,2,2)
                obj.dispV(abs(Vpar),abs(Vparend),'v_{par} [m/s]',[1,timesteppart],'None')
                
                sgtitle(sprintf('t=%1.2e[ns] r=[%2.1f, %2.1f] [mm] z=[%2.1f, %2.1f] [mm]',obj.tpart(timestep)*1e9, rspan*1e3, zspan*1e3))
                obj.savegraph(f,sprintf('%s/%sParts_V_parper',obj.folder,obj.name));
                if gcs
                    obj.savegraph(f,sprintf('%s/%sParts_V_parpergcs',obj.folder,obj.name));
                else
                    obj.savegraph(f,sprintf('%s/%sParts_V_parper',obj.folder,obj.name));
                end
                
            end
            
        end
        
        function display2DVdistrib(obj,timestep, rpos, zpos, gcs)
            %display2DVdistrib plot the velocity distribution function
            % in m/s for the parallel and perpendicular velocity
            % and for the radial azimuthal velocity
            % as a 2D contour plot the show the velocity phase space distribution
            %extracted from the markers at position window from rpos(1)
            %rpos(end) and zpos(1) to zpos(end)
            % and at time obj.tpart(timestep)
            %rpos and zpos are given as grid indices
            % gcs define if you give the perpendicular velocity in the
            % guiding center frame or in the laboratory frame
            if(obj.R.nt>=2)
                if nargin<2
                    timesteppart=length(obj.tpart);
                else
                    timesteppart=timestep;
                end
                if nargin<3 || isempty(rpos)
                    rspan=[obj.rgrid(1) obj.rgrid(end)];
                else
                    r=[obj.rgrid(1);(obj.rgrid(1:end-1)+obj.rgrid(2:end))*0.5;obj.rgrid(end)];
                    rspan=[r(rpos(1)) r(rpos(end)+1)];
                end
                if nargin<4 || isempty(zpos)
                    zspan=[obj.zgrid(1) obj.zgrid(end)];
                else
                    z=[obj.zgrid(1);(obj.zgrid(1:end-1)+obj.zgrid(2:end))*0.5;obj.zgrid(end)];
                    zspan=[z(zpos(1)) z(zpos(end)+1)];
                end
                if nargin<5
                    gcs=false; % define if we look in the guiding center system
                end
                
                nbp=min(obj.nbparts(timesteppart),obj.R.nparts);
                R=obj.R(1:nbp,timesteppart,false);
                Z=obj.Z(1:nbp,timesteppart,false);
                ids=R>=rspan(1) & R<=rspan(2) & Z>=zspan(1) & Z<=zspan(2);
                
                Vperp=obj.Vperp(1:nbp,timesteppart,false,gcs);
                Vpar=obj.Vpar(1:nbp,timesteppart,false);
                Vr=obj.VR(1:nbp,timesteppart,false);
                Vthet=obj.VTHET(1:nbp,timesteppart,false);
                Vper=Vperp(ids);
                Vpar=Vpar(ids);
                
                Vr=Vr(ids);
                Vthet=Vthet(ids);
                nbp=sum(ids(:));
                
                f=figure('Name',sprintf("%s v parper distrib",obj.file));
                subplot(2,1,1)
                [N,Xedges,Yedges] = histcounts2(Vpar,Vper,20);
                Xedges=(Xedges(1:end-1)+Xedges(2:end))/2;
                Yedges=(Yedges(1:end-1)+Yedges(2:end))/2;
                contourf(Xedges,Yedges,N')
                xlabel('v_{par} [m/s]')
                ylabel('v_{\perp} [m/s]')
                c=colorbar;
                c.Label.String='Counts';
                
                subplot(2,1,2)
                [N,Xedges,Yedges] = histcounts2(Vthet,Vr,20);
                Xedges=(Xedges(1:end-1)+Xedges(2:end))/2;
                Yedges=(Yedges(1:end-1)+Yedges(2:end))/2;
                contourf(Xedges,Yedges,N')
                %histogram2(Vthet,Vr,'displaystyle','tile','binmethod','auto')
                %scatter(Vthet,Vr)
                xlabel('v_\theta [m/s]')
                ylabel('v_r [m/s]')
                c=colorbar;
                c.Label.String='Counts';
                
                sgtitle(sprintf('t=%1.2e[ns] r=[%2.1f, %2.1f] [mm] z=[%2.1f, %2.1f] [mm] N=%3i',mean(obj.tpart(timestep))*1e9, rspan*1e3, zspan*1e3,nbp))
                mkdir(sprintf('%s/vdist',obj.folder))
                if gcs
                    obj.savegraph(f,sprintf('%s/vdist/%sParts_V_2dparpergcs_r%iz%it%i',obj.folder,obj.name,floor(mean(rpos)),floor(mean(zpos)),floor(mean(timestep))));
                else
                    obj.savegraph(f,sprintf('%s/vdist/%sParts_V_2dparper_r%iz%it%i',obj.folder,obj.name,floor(mean(rpos)),floor(mean(zpos)),floor(mean(timestep))));
                end
                
            end
            
        end
        
        function [p, maxnb, c]=displayPhaseSpace(obj,type,partsstep, Rindex, Zindex,legendtext, figtitle, f, maxnb, c, gcs)
            if nargin<8
                f=figure;
                f=gca;
            end
            if nargin<7
                figtitle=sprintf('r=%1.2f [mm] z=%1.2f [mm] \\Delta\\phi=%1.1f[kV] R=%1.1f',obj.rgrid(Rindex)*1e3,obj.zgrid(Zindex)*1e3,(obj.potout-obj.potinn)*obj.phinorm/1e3,obj.Rcurv);
            end
            if nargin <6
                legendtext=sprintf('t=%1.3g [s]',obj.tpart(partsstep));
            end
            fieldstep=find(obj.tpart(partsstep(end))==obj.t2d,1);
            
            if nargin>=10
                ctemp=c;
                n=zeros(length(c{1}),length(c{2}));
            else
                nbins=15;
                n=zeros(nbins);
            end
            if nargin <11
                gcs=true;
            end
            
            
            for i=1:length(partsstep)
                odstep=find(obj.tpart(partsstep(i))==obj.t0d);
                nbp=min(obj.R.nparts,obj.nbparts(partsstep(i)));
                Rp=obj.R(1:nbp,partsstep(i),false);
                Zp=obj.Z(1:nbp,partsstep(i),false);
                deltar=obj.dr(2)/2;
                deltarm=obj.rgrid(Rindex)-sqrt(obj.rgrid(Rindex)^2-deltar^2-2*obj.rgrid(Rindex)*deltar);
                deltaz=obj.dz/2;
                Indices=Rp>=obj.rgrid(Rindex)-deltarm & Rp<obj.rgrid(Rindex)+deltar & Zp>=obj.zgrid(Zindex)-deltaz & Zp<obj.zgrid(Zindex)+deltaz;
                
                if strcmp(type,'parper')
                    Vperp=obj.Vperp(Indices,partsstep(i),false,gcs);
                    Vpar=obj.Vpar(Indices,partsstep(i),false,gcs);
                    if exist('ctemp','var')
                        [ntemp,ctemp] = hist3([Vpar, Vperp],'Ctrs',ctemp );
                    else
                        [ntemp,ctemp]=hist3([Vpar, Vperp],nbins*[1 1]);
                    end
                else
                    Vr=obj.VR(Indices,partsstep(i),false);
                    Vthet=obj.VTHET(Indices,partsstep(i),false);
                    if exist('ctemp','var')
                        [ntemp,ctemp] = hist3([Vr, Vthet],'Ctrs',ctemp );
                    else
                        [ntemp,ctemp]=hist3([Vr, Vthet],nbins*[1 1]);
                    end
                end
                
                
                n=(n+ntemp)/2;
            end
            
            c=ctemp;
            
            if nargin<8 || maxnb < 0
                maxnb=max(max(n));
            end
            sumn=sum(n(:))
            contourf(f,c{1},c{2},n'/maxnb,'Displayname',legendtext);
            colorbar(f)
            hold(f,'on')
            
            
            %caxis(f,[0 1])
            
            axis(f, 'equal')
            if strcmp(type,'parper')
                xlabel(f,'v_{par} [m/s]')
                if gcs
                    ylabel(f,'v_\perp^* [m/s]')
                else
                    ylabel(f,'v_\perp [m/s]')
                end
            else
                xlabel(f,'v_r [m/s]')
                ylabel(f,'v_\theta [m/s]')
            end
            xlimits=xlim(f);
            ylimits=ylim(f);
            xtickformat(f,'%.2g')
            ytickformat(f,'%.2g')
            
            
            if strcmp(type,'parper')
                b=obj.rgrid(end);
                a=obj.rgrid(1);
                % drift velocity in vacuum vessel with bias
                vd1=((obj.potout-obj.potinn)*obj.phinorm/obj.rgrid(Rindex))/obj.Bz(Zindex,Rindex)/log(b/a);
                % initial particule velocity ( thermal velocity)
                vinit=sqrt(obj.kb*22000/obj.me);
                
                [Zmesh,~]=meshgrid(obj.zgrid,obj.rAthet(:,Zindex));
                Zeval=obj.zgrid(1:end);
                zleftlim=1;
                zrightlim=length(Zeval);
                Psieval=obj.rAthet(Rindex,Zindex)*ones(length(Zeval),1);
                % Electrostatic potential on magnetic field line
                % coordinates
                phis=griddata(Zmesh,obj.rAthet,obj.pot(:,:,fieldstep),Zeval,Psieval,'natural');
                
                % Magnetic field mirror ratio at each grid position
                % compared to local position
                R=griddata(Zmesh,obj.rAthet,obj.B',Zeval,Psieval,'natural')/obj.B(Zindex,Rindex);
                
                
                if ~gcs
                    [~,tfield]=min(abs(obj.t2d-obj.tpart(partsstep(1))));
                    timeEr=obj.Er(:,:,tfield);
                    timeEz=obj.Ez(:,:,tfield);
                    posindE=sub2ind(size(timeEr),Rindex,Zindex);
                    % ExB drift velocity at considered time step
                    Vdrift=(timeEz(posindE).*obj.Br(Zindex,Rindex)-timeEr(posindE).*obj.Bz(Zindex,Rindex))./obj.B(Zindex,Rindex).^2;
                else
                    Vdrift=0;
                end
                Rcurvl=R(zleftlim);
                Rcurvr=R(zrightlim);
                deltaphicompl=0.5*obj.me/obj.qe*((vd1+vinit)^2*(Rcurvl-1));
                deltaphil=-phis(Zindex)+phis(zleftlim);
                deltaphicompr=0.5*obj.me/obj.qe*((vd1+vinit)^2*(Rcurvr-1));
                deltaphir=-phis(Zindex)+phis(zrightlim);
                
                vpar=linspace(1,3*max(abs(Vpar)),1000);
                vperstarr=sqrt((2*obj.qe/obj.me*deltaphir+vpar.^2)/(Rcurvr-1));
                vperr=sqrt(vperstarr.^2+Vdrift^2-2*Vdrift*abs(vperstarr));
                vparr=vpar(real(vperr)~=0& imag(vperr)==0);
                vperr=vperr(real(vperr)~=0& imag(vperr)==0);
                
                vperstarl=sqrt((2*obj.qe/obj.me*deltaphil+vpar.^2)/(Rcurvl-1));
                vperl=sqrt(vperstarl.^2+Vdrift^2-2*Vdrift*abs(vperstarl));
                vparl=vpar(real(vperl)~=0 & imag(vperl)==0);
                vperl=vperl(real(vperl)~=0 & imag(vperl)==0);
                p=plot(f,vparr,vperr,'r-','displayname','Simulation','linewidth',2);
                plot(f,-vparl,vperl,'r--','linewidth',2)
                
                % Prediction using model
                vpar=linspace(1,3*max(abs(Vpar)),1000);
                vperstarr=sqrt((2*obj.qe/obj.me*deltaphicompr+vpar.^2)/(Rcurvr-1));
                vperr=sqrt(vperstarr.^2+Vdrift^2-2*Vdrift*abs(vperstarr));
                vparr=vpar(real(vperr)~=0);
                vperr=vperr(real(vperr)~=0);
                
                vperstarl=sqrt((2*obj.qe/obj.me*deltaphicompl+vpar.^2)/(Rcurvl-1));
                vperl=sqrt(vperstarl.^2+Vdrift^2-2*Vdrift*abs(vperstarl));
                vparl=vpar(real(vperl)~=0);
                vperl=vperl(real(vperl)~=0);
                p=plot(f,vparr,vperr,'r-.','displayname','Prediction','linewidth',2);
                plot(f,-vparl,vperl,'r:','linewidth',2)
                
            else
                p=gobjects(0);
            end
            
            xlim(f,xlimits)
            ylim(f,ylimits)
            legend(p(isgraphics(p)))
            legend(p(isgraphics(p)),'location','north')
            
            title(f,figtitle)
            
        end
        
        function displaytavgdensity(obj,fieldstart,fieldend)
            %displaytavgdensity plot the fluid density averaged in time for times t2d(fieldstart:fieldend)
            % Displays also the geometry of the metallic boundaries and the
            % magnetic field lines
            f=figure('Name', sprintf('%sfluid_dens',obj.name));
            ax1=gca;
            geomw=obj.geomweight(:,:,1);
            dens=mean(obj.N(:,:,fieldstart:fieldend),3);
            dens(geomw<=0)=0;
            maxdens=max(dens(:));
            [h,curve]=contourf(ax1,obj.zgrid*1000,obj.rgrid*1000,dens,30);
            hold on
            set(curve,'linecolor','none');
            %contour(ax1,obj.zgrid,obj.rgrid,obj.geomweight(:,:,1),[0 0],'r-','linewidth',1.5,'Displayname','Boundaries');
            if(obj.conformgeom)
                ylim(ax1,[obj.rgrid(1) obj.rgrid(sum(obj.nnr(1:2)))]*1000)
            else
                ylim(ax1,[obj.rgrid(1) obj.rgrid(end)]*1000)
            end
            xlim(ax1,[obj.zgrid(1) obj.zgrid(end)]*1000)
            xlabel(ax1,'z [mm]')
            ylabel(ax1,'r [mm]')
            %title(ax1,sprintf('Density t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',obj.t2d(fieldstart),obj.t2d(fieldend),double(maxdens)))
            c = colorbar(ax1);
            c.Label.String= 'electron density [m^{-3}]';
            view(ax1,2)
            colormap(flipud(hot));
            
            % Draws the magnetic field lines
            Blines=obj.rAthet;
            levels=linspace(min(Blines(obj.geomweight(:,:,1)>0)),max(Blines(obj.geomweight(:,:,1)>0)),20);
            Blines(obj.geomweight(:,:,1)<0)=NaN;
            [~,h1]=contour(ax1,obj.zgrid*1000,obj.rgrid*1000,Blines,real(levels),'-.','color','k','linewidth',1.5,'Displayname','Magnetic field lines');
            
            % Draw the metallic boundaries and the geometry itself
            [c1,hContour]=contourf(ax1,obj.zgrid*1000,obj.rgrid*1000,-geomw,[0,0],'linewidth',1.5);
            drawnow;
            xlim(ax1,[obj.zgrid(1)*1000 obj.zgrid(end)*1000])
            % Change the color of the metallic boundaries to grey
            hFills=hContour.FacePrims;
            [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
            try
                hFills(end).ColorData = uint8([150;150;150;255]);
                for idx = 1 : numel(hFills)-1
                    hFills(idx).ColorData(4) = 0;   % default=255
                end
            catch
            end

            grid on;
            hold on;
            f.PaperOrientation='landscape';
            f.PaperUnits='centimeters';
            papsize=[16 14];
            f.PaperSize=papsize;
            set(ax1,'fontsize',14)
            %axis equal
            obj.savegraph(f,sprintf('%sfluid_dens',obj.name))
        end
        
        function displayconfiguration(obj,fieldstart,fieldend)
            %displayconfiguration plot the configuration of the simulation
            % domain withe boundaries, the magnetic field lines the
            % electric equipotential lines and the electron density
            % averaged in time between t2d(fieldstart) and t2d(fieldend)
            dens=mean(obj.N(:,:,fieldstart:fieldend),3);
            geomw=obj.geomweight(:,:,1);
            maxdens=max(dens(:));
            geomw(geomw<0)=0;
            geomw(geomw>0)=maxdens;
            dens(geomw<=0)=0;
            geomw(geomw>0)=NaN;
            
            
            
            f=figure('Name', sprintf('%s fields',obj.name));
            ax1=gca;
            title(sprintf('Configuration'))
            %dens(dens<=1e13)=NaN;
            %% electron density
            h=contourf(ax1,obj.zgrid*1000,obj.rgrid*1000,dens,50,'Displayname','n_e [m^{-3}]', 'linestyle','none');
            hold on;
            colormap(flipud(hot));
            %% Magnetic field lines
            Blines=obj.rAthet;
            levels=linspace(min(Blines(obj.geomweight(:,:,1)>0)),max(Blines(obj.geomweight(:,:,1)>0)),20);
            [~,h1]=contour(ax1,obj.zgrid*1000,obj.rgrid*1000,Blines,real(levels),'-.','color','k','linewidth',1.5,'Displayname','Magnetic field lines');
            
            %% Equipotential lines
            Pot=mean(obj.pot(:,:,fieldstart:fieldend),3);
            Pot(obj.geomweight(:,:,1)<0)=NaN;
            %levels=8;%[-3.4 -5 -10 -15 -20 -25];%7;
            potcolor='b';%[0.3660 0.6740 0.1880];
            [c1,h2]=contour(ax1,obj.zgrid*1000,obj.rgrid*1000,Pot,'--','color',potcolor,'ShowText','on','linewidth',1.2,'Displayname','Equipotentials [kV]');
            clabel(c1,h2,'Color',potcolor)
            
            % Grey outline shows metallic parts
            [c1,hContour]=contourf(ax1,obj.zgrid*1000,obj.rgrid*1000,geomw, [0 0]);
            
            drawnow;
            
            % set the axia limits
            xlim(ax1,[obj.zgrid(1)*1000 obj.zgrid(end)*1000])
            if(obj.conformgeom)
                ylim([ax1 ],[obj.rgrid(1)*1000 obj.rgrid(rgridend)*1000])
            else
                ylim([ax1],[obj.rgrid(1)*1000 obj.rgrid(end)*1000])
            end
            legend([h1,h2],{'Magnetic field lines','Equipotentials [V]'},'location','northeast','Interpreter','latex')
            xlabel(ax1,'Z [mm]','Interpreter','latex')
            ylabel(ax1,'R [mm]','Interpreter','latex')
            set(legend,'FontSize',18);
            set (gca, 'fontsize', 22)

            
            c = colorbar(ax1);
            c.Label.String= 'n[m^{-3}]';
            view(ax1,2)
            
            grid on;
            hFills=hContour.FacePrims;
            [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
            try
                hFills(1).ColorData = uint8([150;150;150;255]);
                for idx = 2 : numel(hFills)
                    hFills(idx).ColorData(4) = 0;   % default=255
                end
            catch
            end
            [~, name, ~] = fileparts(obj.file);
            
            % with this you could show the outline of the maxwellian source
            % if obj.maxwellsrce.present
            %     rlen=diff(obj.maxwellsrce.rlim);
            %     zlen=diff(obj.maxwellsrce.zlim);
            % rectangle('Position',[obj.maxwellsrce.zlim(1) obj.maxwellsrce.rlim(1) zlen rlen]*1000,'Edgecolor','g','Linewidth',2,'Linestyle','--')
            % end
            
            % in case of coaxial configuration, extend the display domain
            % and add grey rectangles to show metallic boundaries
            if( obj.walltype >=2 && obj.walltype<=4)
                rectangle('Position',[obj.zgrid(1) obj.r_b obj.zgrid(end)-obj.zgrid(1) 0.001]*1e3,'FaceColor',[150 150 150]/255,'Edgecolor','none')
                ylimits=ylim;
                ylim([ylimits(1),ylimits(2)+1])
            end
            if sum(obj.geomweight(:,1,1))==0
                rectangle('Position',[obj.zgrid(1) obj.r_a-0.001 obj.zgrid(end)-obj.zgrid(1) 0.001]*1e3,'FaceColor',[150 150 150]/255,'Edgecolor','none')
                ylimits=ylim;
                ylim([ylimits(1)-1,ylimits(2)])
            end
            f.PaperUnits='centimeters';
            %axis equal
            
            papsize=[14 5 ];
            
            
            obj.savegraph(f,sprintf('%s/%sFields',obj.folder,obj.name),papsize);
        end
        
        function displaymagfield(obj)
            %displaymagfield display the magnetic field lines and the
            %amplitude of the magnetic field using a contour
            % also show the domain boundaries
            
            B=obj.B';
            
            f=figure('Name', sprintf('%s B field',obj.name));
            B(obj.geomweight(:,:,1)<0)=NaN;
            ax1=gca;
            title(sprintf('Configuration'))
            h=contourf(ax1,obj.zgrid*1000,obj.rgrid*1000,B,linspace(min(B(:)),max(B(:)),50),'Displayname','B [T]', 'linestyle','none');
            hold on;
            
            %% Magnetic field lines
            Blines=obj.rAthet;
            levels=linspace(min(Blines(obj.geomweight(:,:,1)>0)),max(Blines(obj.geomweight(:,:,1)>0)),50);
            Blines(obj.geomweight(:,:,1)<0)=NaN;
            [~,h1]=contour(ax1,obj.zgrid*1000,obj.rgrid*1000,Blines,real(levels),'r-.','linewidth',1.5,'Displayname','Magnetic field lines');
            
            colormap(ax1,'parula')
            
            % Grey outline
            geomw=obj.geomweight(:,:,1);
            geomw(geomw>0)=NaN;
            geomw(geomw<0)=min(B(:));
            [c1,hContour]=contourf(ax1,obj.zgrid*1000,obj.rgrid*1000,geomw, [0 0]);
            
            drawnow;
            xlim(ax1,[obj.zgrid(1)*1000 obj.zgrid(end)*1000])
            if(obj.conformgeom)
                ylim([ax1 ],[obj.rgrid(1)*1000 obj.rgrid(rgridend)*1000])
            else
                ylim([ax1],[obj.rgrid(1)*1000 obj.rgrid(end)*1000])
            end
            legend([h1],{'Magnetic field lines'},'location','northwest')
            xlabel(ax1,'z [mm]')
            ylabel(ax1,'r [mm]')
            %title(ax1,sprintf('Density t=[%1.2g-%1.2g]s n_e=%1.2gm^{-3}',M.t2d(fieldstart),M.t2d(fieldend),double(maxdens)))
            c = colorbar(ax1);
            c.Label.String= 'B [T]';
            view(ax1,2)
            %set(h,'edgecolor','none');
            grid on;
            hFills=hContour.FacePrims;
            [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
            %caxis([min(B(:)) max(B(:))])
            try
                hFills(1).ColorData = uint8([150;150;150;255]);
                for idx = 2 : numel(hFills)
                    hFills(idx).ColorData(4) = 0;   % default=255
                end
            catch
            end
            [~, name, ~] = fileparts(obj.file);
            if( obj.walltype >=2 && obj.walltype<=4)
                rectangle('Position',[obj.zgrid(1) obj.r_b obj.zgrid(end)-obj.zgrid(1) 0.001]*1e3,'FaceColor',[150 150 150]/255,'Edgecolor','none')
                ylimits=ylim;
                ylim([ylimits(1),ylimits(2)+1])
            end
            if(isempty(obj.spl_bound))
                rectangle('Position',[obj.zgrid(1) obj.r_a-0.001 obj.zgrid(end)-obj.zgrid(1) 0.001]*1e3,'FaceColor',[150 150 150]/255,'Edgecolor','none')
                ylimits=ylim;
                ylim([ylimits(1)-1,ylimits(2)])
            end
            f.PaperUnits='centimeters';
            %axis equal
            
            papsize=[14 5 ];
            pos=f.Position;
            pos(3)=floor(1.3*pos(3));
            f.Position=pos;
            
            obj.savegraph(f,sprintf('%s/%s_Bfield',obj.folder,obj.name),papsize);
        end
        
        function displaySurfFlux(obj,timestep, subdiv)
            %displaySurfFlux plot the current densities
            %on the domain boundaries for time t2d(timestep)
            %directly on the boundaries themselves
            %make it easier to see where the currents are collected
            
            if nargin<3
                subdiv=1;
            end
            mflux= obj.Metallicflux(timestep,subdiv);
            lflux= -squeeze(obj.Axialflux(timestep,1))';
            rflux= squeeze(obj.Axialflux(timestep,length(obj.zgrid)))';
            
            time=obj.t2d(timestep);
            
            if nargin<3
                ids=1:length(mflux);
            end
            
            %%
            P=obj.neutcol.neutdens*obj.kb*300/100;% pressure at room temperature in mbar
            f=figure('name','fluxevol');
            linew=5;
            %obj.displaysplbound(gca,1e3);
            contour(obj.zgrid*1e3,obj.rgrid*1e3,obj.geomweight(:,:,1),[0 0],'b-','linewidth',1.5);
            hold on
            for i=1:length(mflux.p)
                x=mflux.p{i}(1,:)*1000;
                y=mflux.p{i}(2,:)*1000;
                y(end)=NaN;
                c=mflux.gamma{i}'*obj.qe/(100^2)/P;
                c(c<=0)=NaN;
                patch(x,y,c,'EdgeColor','interp','LineWidth',linew);
                hold on
            end
            
            x=obj.zgrid(1)*ones(size(obj.rgrid))*1000;
            y=obj.rgrid*1000;
            y(end)=NaN;
            c=lflux*obj.qe/(100^2)/P;
            c(c<=0)=NaN;
            patch(x,y,c,'EdgeColor','interp','LineWidth',linew);
            
            
            x=obj.zgrid(end)*ones(size(obj.rgrid))*1e3;
            y=obj.rgrid*1000;
            y(end)=NaN;
            c=rflux*obj.qe/(100^2)/P;
            c(c<=0)=NaN;
            patch(x,y,c,'EdgeColor','interp','LineWidth',linew);
            
            title(sprintf('t=%4.2f [ns]',time*1e9))
            
            c=colorbar;
            c.Label.String= 'j\cdotn [A/(cm^2 mbar)]';
            xlabel('z [mm]')
            ylabel('r [mm]')
            colormap(jet)
            set(gca,'colorscale','log')
            
            %% Magnetic field lines
            Blines=obj.rAthet;
            levels=linspace(min(Blines(obj.geomweight(:,:,1)>0)),max(Blines(obj.geomweight(:,:,1)>0)),20);
            Blines(obj.geomweight(:,:,1)<0)=NaN;
            [~,h1]=contour(obj.zgrid*1000,obj.rgrid*1000,Blines,real(levels),'m-.','linewidth',1.5,'Displayname','Magnetic field lines');
            %axis equal
            
            obj.savegraph(f,sprintf('%s/%s_surfFlux_it2d_%i',obj.folder,obj.name,timestep),[16 14]);
        end
        
        % interactive window to display the terms of the pressure tensor
        dispespicPressure(obj,logdensity,showgrid,fixed,temperature)
        
        %  interactive window to display the electron density, magnetic
        %  field lines, electric potential and field at given time steps
        dispespicFields(obj,logdensity,showgrid,fixed,parper)
        
        
        function displaycollfreq(obj)
            %displaycollfreq plot the collision frequencies in Hz/mbar for a range of
            %electron kinetic energies in eV for the different collision
            %processes considered: ionisation and elastic collisions
            E=logspace(1,4,1000);
            
            v=sqrt(2/obj.msim*obj.weight*E*obj.qe);
            P=obj.neutcol.neutdens*obj.kb*300/100;% pressure at room temperature in mbar
            tauio=P./(obj.neutcol.neutdens*obj.sigio(E).*v);
            tauiom=P./(obj.neutcol.neutdens*obj.sigmio(E).*v);
            
            tauelam=P./(obj.neutcol.neutdens*obj.sigmela(E).*v);
            
            f=figure('name','t scales coll');
            loglog(E,1./tauio,'displayname','ionisation','linewidth',1.5)
            hold on
            loglog(E,1./tauiom,'displayname','ionisation momentum','linewidth',1.5)
            loglog(E,1./tauelam,'displayname','elastic','linewidth',1.5)
            loglog(E,1./tauio+1./tauiom+1./tauelam,'displayname','total drag','linewidth',1.5)
            
            loglog([E(1) E(end)],1./(2*pi/obj.omece)* [1 1],'--','displayname','cyclotronic','linewidth',1.5)
            
            xlabel('Electron kinetic energy [eV]')
            ylabel('\nu [Hz/mbar]')
            legend('location','southeast')
            grid on
            
            obj.savegraph(f,sprintf('%s/collfreqscales',obj.folder),[14 12]);
        end
        
        function displaycrosssec(obj)
            %displaycrosssec plot the collision crosssections in m^2 for a range of
            %electron kinetic energies in eV for the different collision
            %processes considered: ionisation and elastic collisions
            E=logspace(1,4,1000);
            
            
            sig_io=obj.sigio(E);
            sig_iom=obj.sigmio(E);
            
            sig_elam=obj.sigmela(E);
            
            f=figure('name','t scales coll');
            loglog(E,sig_io,'displayname','ionisation','linewidth',1.5)
            hold on
            loglog(E,sig_iom,'displayname','ionisation momentum','linewidth',1.5)
            loglog(E,sig_elam,'displayname','elastic','linewidth',1.5)
            loglog(E,sig_io+sig_elam+sig_iom,'displayname','total drag','linewidth',1.5)
            
            
            xlabel('Energy [eV]')
            ylabel('\sigma [m^{2}]')
            legend('location','southeast')
            grid on
            
            obj.savegraph(f,sprintf('%s/coll_cross_sec_scales',obj.folder),[14 12]);
        end
        
        %------------------------------------------
        %  Helper functions needed for other functions
        
        function [zpos,rpos]=getpos(obj,tstep)
            % interactive window to return an specific axial and radial
            % position picked from the cloud density
            if nargin<2
                tstep=length(obj.t2d);
            end
            n=obj.N(:,:,tstep);
            n(obj.geomweight(:,:,1)<0)=NaN;
            
            figure
            contourf(obj.zgrid,obj.rgrid,n);
            xlabel('z [m]')
            ylabel('r [m]')
            [x,y]=ginput(1);
            zpos=find(x>obj.zgrid,1,'last');
            rpos=find(y>obj.rgrid,1,'last');
            
            hold on
            plot(obj.zgrid(zpos),obj.rgrid(rpos),'rx')
            
            fprintf('zpos=%i  rpos=%i  z=%1.4f r=%1.4f\n',zpos,rpos,obj.zgrid(zpos),obj.rgrid(rpos))
            
        end
        
        function changed=ischanged(obj)
            %ischanged Check if the file has been changed since the initial loading of the file
            %and if some data must be reloaded
            try
                filedata=dir(obj.fullpath);
                checkedtimestamp=filedata.date;
                if (max(checkedtimestamp > obj.timestamp) )
                    changed=true;
                    return
                end
                changed=false;
                return
            catch
                changed=true;
                return
            end
        end
        
        function dispV(obj,V,Vend,label,t, dist, vd)
            %dispV generic functio to plot the velocity distribution and
            %comapare two timesteps V and Vend at time  t(1) and t(2)
            if nargin<6
                dist='gaussian';
            end
            if nargin<7
                vd=0;
            end
            vmean=mean(V(~isnan(V)));
            vtherm=std(V(~isnan(V)),1);
            vmeanend=mean(Vend(~isnan(Vend)));
            vthermend=std(Vend(~isnan(Vend)),1);
            
            if(length(V)>1)
                [Counts,edges]=histcounts(V,'binmethod','sqrt');
                binwidth=mean(diff(edges));
                plot([edges(1) 0.5*(edges(2:end)+edges(1:end-1)) edges(end)],[0 Counts 0],'DisplayName',sprintf("t=%2.3d [ns]",obj.tpart(t(1))*1e9));
                hold on
            end
            hold on
            [Counts,edges]=histcounts(Vend,'binmethod','sqrt');
            plot([edges(1) 0.5*(edges(2:end)+edges(1:end-1)) edges(end)],[0 Counts 0],'DisplayName',sprintf("t=%2.3d [ns]",obj.tpart(t(2))*1e9));
            
            if strcmp(dist,'maxwell')
                vfit=linspace(0,edges(end),300);
                a=vmeanend/sqrt(2);
                dist=sqrt(2/pi)*vfit.^2.*exp(-((vfit).^2-vd^2)/2/a^2)/a^3;
                dist=dist/max(dist);
                plot(vfit,max(Counts)*dist,'displayname',sprintf('Maxw mu=%2.2g sigma=%2.2g',vmeanend,vthermend))
            elseif strcmp(dist,'gaussian')
                vfit=linspace(edges(1),edges(end),300);
                dist=exp(-(vfit-vmeanend).^2/2/vthermend^2);
                plot(vfit,max(Counts)*dist,'displayname',sprintf('gauss mu=%2.2g sigma=%2.2g',vmeanend,vthermend))
            end
            ylabel('counts')
            xlabel(label)
            grid on
            legend('location','southoutside','orientation','vertical')
        end
        
        function cross_sec=fit_cross_sec(obj,energy,crosssec_table)
            %Interpolate the cross-section at the given energy using the
            %crosssec_table and an exponential fitting
            cross_sec=0;
            if (energy<=0 || isnan(energy) || isinf(energy))
                return
            end
            id=find(energy>crosssec_table(:,1),1,'last');
            if(isempty(id))
                id=1;
            end
            id=min(size(crosssec_table,1)-1,id);
            id=max(1,id);
            
            cross_sec=crosssec_table(id,2)*(energy/crosssec_table(id,1))^crosssec_table(id,3);
        end
        
        function fighandle=savegraph(obj, fighandle, name, papsize)
            %% Saves the given figure as a pdf a .fig and an eps using export_fig
            fighandle.PaperUnits='centimeters';
            if (nargin < 4)
                papsize=[14 16];
            end
            set(fighandle, 'Color', 'w');
            fighandle.PaperSize=papsize;
            %export_fig(fighandle,name,'-png','-r300')
            exportgraphics(fighandle,sprintf('%s.png',name),'Resolution',300)
            print(fighandle,name,'-dpdf','-fillpage')
            savefig(fighandle,name)
            set(fighandle, 'Color', 'w');
            exportgraphics(fighandle,sprintf('%s.eps',name))
            %export_fig(fighandle,name,'-eps','-painters')
            
        end
        
        function sig=dsigmaio(obj,Ekin, Ebar, Ei, E0, chi, gamma)
            % calculates the integrand used for the ionisation collision
            % cross section for momentum exchange for the incoming electron
            % it is only used by obj.sigmiopre
            gamma=reshape(gamma,1,[],1);
            chi=reshape(chi,1,1,[]);
            
            siggamma=sin(gamma).*(E0^2+8*(1-chi)*(Ekin-Ei)*E0)./(E0+4*(1-chi)*(Ekin-Ei)-4*(1-chi)*(Ekin-Ei).*cos(gamma)).^2/2;
            
            sigchi=(Ekin-Ei)./(Ebar*atan((Ekin-Ei)/(2*Ebar)).*(1+(chi*(Ekin-Ei)/Ebar).^2));
            
            dp=1- trapz(gamma,sqrt((1-chi).*(1-Ei/Ekin)).*cos(gamma).*siggamma,2);%- trapz(gamma,sqrt((1-chi).*(1-Ei/Ekin)).*cos(gamma).*siggamma,2);
            sig=sigchi.*dp;
        end
        
        function sigm=sigmiopre(obj,E, init)
            % returns the precalculated values used for the interpolation
            % of the ionisation collision cross-section for momentum
            % exchange for the incoming electron
            if nargin <3
                init=false;
            end
            if(~init &&( ~obj.neutcol.present || isempty(obj.neutcol.io_cross_sec)))
                sigm=zeros(size(E));
                return
            end
            Ebar=obj.neutcol.scatter_fac;
            Ei=obj.neutcol.Eion;
            E0=obj.neutcol.E0;
            nE=numel(E);
            
            nchi=300;
            ngamma=300;
            gamma=linspace(0,pi,ngamma);
            chi=linspace(0,0.5,nchi);
            %sigm2=zeros(nE,nchi);
            sigm=zeros(size(E));
            
            for i=1:nE
                if(E(i)>=Ei)
                    sigm2=zeros(nchi,1);
                    for j=1:nchi
                        %sigm2(j)=trapz(alpha,trapz(gamma,obj.dsigmaio(E(i),Ebar,Ei,E0,chi(j),alpha,gamma),2),1);
                        sigm2(j)=obj.dsigmaio(E(i),Ebar,Ei,E0,chi(j),gamma);
                    end
                    sigm(i)=trapz(chi,sigm2)*obj.sigio(E(i),init);
                    %sigm(i)=trapz(chi,trapz(alpha,trapz(gamma,dsigmaio(obj,E(i),Ebar,Ei,E0,chi,alpha,gamma),2),1),3)*obj.sigio(E(i),init);
                end
            end
        end
        
    end
end
