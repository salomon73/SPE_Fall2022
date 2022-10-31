%function initialiseMaxwellVzVr(filename, lowerBound, upperBound, nPoints,TextFile)
    %% Define the constants %%
    kB = (25.7/298)*1e-3;      % eV/K
    m  = 9.10938300000000e-31; % electron mass
    e  = 1.60217662000000e-19; % J/eV
    
    %% Create geom structure %% 
    %geom = espic2dhdf5(filename);

    lowerBound = 0.1;
    upperBound = 20;
    nPoints = 20;
    cathode = geom.spl_bound.boundary(1);
    LowPts = 15;
    UpPts  = 65;
    
    f   = cathode.fun;
    E   = linspace(lowerBound, upperBound, nPoints);
    vR  = sqrt((2*e/m)*E);
    Points = zeros(UpPts-LowPts,2);
    Tang   = zeros(UpPts-LowPts,2);
    Norm   = zeros(UpPts-LowPts,2);
    df = fnder(f,1);
    
    for ii = 0:UpPts-LowPts-1
        
       pts = spval(f,f.knots(ii+LowPts));
       Points(ii+1,1) = pts(1);
       Points(ii+1,2) = pts(2);
       tgt = spval(df,df.knots(ii+LowPts-1));
       Tang(ii+1,1) =  tgt(1);
       Tang(ii+1,2) =  tgt(2);
       Norm(ii+1,1) = -tgt(2);
       Norm(ii+1,2) =  tgt(1);
    end


    LowInd  = LowPts-1;
    HighInd = UpPts-1;
        figure
        hold on 
        plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
        ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
        xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
        set (gca, 'fontsize', 20)
        hold on 
        plot(Points(1:3:end,1),Points(1:3:end,2),'r*')
        hold on 
        quiver(Points(1:3:end,1),Points(1:3:end,2), Norm(1:3:end,1), Norm(1:3:end,2), 'b')
        
        
        
