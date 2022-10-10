%% To work from home on local machine (macbook) 
cd /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/Inputs/Test_ions
addpath /Users/salomonguinchard/Documents/GitHub/SPE_Fall2022/matlab_routines
Ions = espic2dhdf5('stable_13_fine.h5');

%% 
Ions_mass = 3.347e-27;
Vr = Ions.VR;
Vz = Ions.VZ;
Vt = Ions.VTHET;
VR = Vr(:,:);
VZ = Vz(:,:);
VT = Vt(:,:);

R = Ions.R;
Z = Ions.Z;

Ra = 0.001;
Rb = 0.01;
Za = -0.32;
Z_b = 0.32;

dt = Ions.dt;
nt = R.nt;
t  = 5*dt*linspace(0,200,nt);
dr = 1.75e-3;

npart = R.nparts;
mat = zeros(npart,nt);
R = R(:,:);

for jj = 1:npart
    for ii = 1:nt
        if R(2,ii) < Ra+dr
            disp 'touch';
            mat(jj,ii) = ii; % Trouve l'indice des particles et du temps auquel elles se trouvent pres de/sur l'électrode 
        else 
            mat(jj,ii) = 0;
            
        end
    end
end

indices = zeros(1,2);
for jj =1:npart
    for ii = 1:nt
        t_ind = find(mat(jj,:)); % trouve le premier indice de temps où il y'a contact pour la particule jj
        indices(jj,1) = jj;
        indices(jj,2) = t_ind(1);
    end
end
%%
k=150; % Index of particle of interest
EnergyIon_n = (1/1.602e-19)*0.5*Ions_mass*(VR(k,indices(k,2)).^2+VZ(k,indices(k,2)).^2+VT(k,indices(k,2)).^2);

for n=1:npart
    E(n) = (1/1.602e-19)*0.5*Ions_mass*(VR(n,indices(n,2)).^2+VZ(n,indices(n,2)).^2+VT(n,indices(n,2)).^2);
end

figure
plot(linspace(1,1000,1000),E, 'k+')
%% Display particles data
dispespicParts(Ions);