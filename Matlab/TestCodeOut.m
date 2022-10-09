%% To work from home on local machine (macbook) 
cd /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Inputs/Test_ions
addpath /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/matlab_routines
Ions = espic2dhdf5('stable_13_fine.h5');

%% 

Vr = Ions.VR;
Vz = Ions.VZ;

R = Ions.R;
Z = Ions.Z;

Ra = 0.001;
Rb = 0.01;
Za = -0.32;
Z_b = 0.32;

dt = Ions.dt;
nt = R.nt;
t  = 5*dt*linspace(0,200,nt);
dr = 1e-3;

npart = R.nparts;
mat = zeros(npart,nt);

for jj = 1:npart
    for ii = 1:nt
        if R(jj,ii) > Rb-dr
            disp 'touch';
            mat(jj,ii) = ii;
        else 
            mat(jj,ii) = NaN;
            
        end
    end
end

%% Display particles data
dispespicParts(Ions);