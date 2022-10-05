

E1 = [.02 .03 .04 .05 .06 .07 .08 .09 ];
E2 = [.1 .2 .3 .4 .5 .6 .7 .8 .9 ];
E3 = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5];
E4 = [10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95];
E5 = [100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 10000];



E = cat(2, E1, E2, E3, E4, E5);

Eloss1 = [958.52 1103.4 1177.2 1209.7 1218.5 1213.4 1200.5 1183.1];
Eloss2 = [1163.3 969.92 838.72 747.47 678.39 622.57 575.31 530.43 493.76];
Eloss3 = [463.10 359.80 297.61 255.35 224.62 201.17 182.63 167.55 155.02 ...
           144.44 135.40 127.52 120.60 114.46 108.99 104.07 99.622 95.579];
Eloss4 = [91.884 67.098 53.524 44.910 38.903 34.461 31.033 28.303 26.075 ...
           24.219 22.649 21.301 20.132 19.106 18.202 17.396 16.673 16.021 ];
Eloss5 = [15.430 11.584 9.5581 8.3666 7.5440 6.9545 6.5131 6.1718 5.9011 ... 
           5.6822 5.5024 5.3527 5.2267 5.1197 5.0281 4.9492 4.8809 4.8124 4.7694];

Eloss = cat(2, Eloss1, Eloss2, Eloss3, Eloss4, Eloss5);


%% 
figure
semilogx(E,Eloss, 'r+-', 'linewidth', 2 )
ylabel('$\frac{dE}{dx}$ [Mev/cm]', 'interpreter', 'latex','Fontsize', 18)
xlabel('$E$ [Mev]', 'interpreter', 'latex', 'Fontsize', 18)


%%
cd /home/sguincha/espic2d/wk/Test_ions
addpath '/home/sguincha/espic2d/matlab/'
Ions = espic2dhdf5('stable_13_fine.h5');

%% 
dispespicParts(Ions)

%% Compute ion energy

Vr = Ions.VR;
Vz = Ions.VZ;
Vt = Ions.VTHET;
Ions_mass = 3.347e-27;
dt = Ions.dt;
t = dt*linspace(0,1000,1000);
EnergyIon1 = 0.5*Ions_mass*(Vr(1:1000,1).^2+Vz(1:1000,1).^2+Vt(1:1000,1).^2);


figure
plot(t,EnergyIon1,'k+-')
xlabel('$t$ [s]', 'interpreter', 'latex','Fontsize', 18)
ylabel('$E$ [eV]', 'interpreter', 'latex', 'Fontsize', 18)
