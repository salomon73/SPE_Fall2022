%function initialiseMaxwellVzVr(filename, lowerBound, upperBound, nPoints,TextFile)
    %% Define the constants %%
    kB = (25.7/298)*1e-3;      % eV/K
    m  = 9.10938300000000e-31; % electron mass
    e  = 1.60217662000000e-19; % J/eV
    
    %% Create geom structure %% 
    %geom = espic2dhdf5('result_46dA_25kv_10mubar.h5');

    lowerBound = 0.1;
    upperBound = 20;
    nPoints = 10;
    cathode = geom.spl_bound.boundary(1);
    LowPts  = 15;
    UpPts   = 65;
    
    f   = cathode.fun;
    E   = linspace(lowerBound, upperBound, nPoints);
    Points = zeros(UpPts-LowPts,2);
    Tang   = zeros(UpPts-LowPts,2);
    Norm   = zeros(UpPts-LowPts,2);
    df = fnder(f,1);
    textFile = 'Test';
    
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

    for jj =1:length(Norm(:,1))

        Norm(jj,:) = 1/norm(Norm(jj,:),2)*Norm(jj,:);   
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
        axis equal
        
        
    fileId = fopen(strcat(textFile,'.txt'),'w');
    fprintf(fileId,'//parts\n');
    for ii =1:length(E)
        for jj =1:UpPts-LowPts
            fprintf(fileId,'%.8f %.1f %.8f %.8f %.1f %.8f \n', Points(jj,1), 0.0, Points(jj,2), Norm(jj,1), 0.0 ,Norm(jj,2));
        end
    end
    fclose(fileId);
    
    
    %%

format long
nppts = UpPts-LowPts;
vR    = zeros(length(E),nppts); 
vZ    = zeros(length(E),nppts); 

for kk=1:length(E)    
    vRmax = sqrt(2*E(kk)/m*e);
    vR(kk,:)    = linspace(0,vRmax, nppts);
    
    for j = 1 :length(vR(kk,:))
        if 2*E(kk)/m*e - vR(kk,j)^2 < 1e-3
            vZ(kk,j) = 0.0;
        else 
            vZ(kk,j) = sqrt(2*E(kk)/m*e - vR(kk,j)^2);
        end
    end
end    

V0R = zeros(length(E),nppts);
V0Z = zeros(length(E),nppts);

for ii =1:length(E)
    for jj=1:nppts
        V0R(ii,jj) = Norm(jj,2)*vR(ii,jj);
        V0Z(ii,jj) = Norm(jj,1)*vZ(ii,jj);
    end
end

    figure
        hold on 
        plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
        ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
        xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
        set (gca, 'fontsize', 20)
        hold on 
        plot(Points(1:3:end,1),Points(1:3:end,2),'r*')
        hold on 
        quiver(Points(1:3:end,1),Points(1:3:end,2), V0Z(1,1:3:end)', V0R(1,1:3:end)', 'b')
        axis equal

    figure
        plot(1:nppts,1/(e*2)*m*(vR(1,:).^2+vZ(1,:).^2),'linewidth', 2)
        hold on 
        plot(1:nppts,1e-5*vR(1,:), 'linewidth', 2)
        hold on 
        plot(1:nppts,1e-5*vZ(1,:), 'linewidth', 2)
        legend('E', 'vR','vZ', 'Location','northwest','Interpreter','latex');
        set(legend,'FontSize',18);
        set (gca, 'fontsize', 22)

%     figure
%         subplot(2,2,1)
%             hold on 
%             plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
%             ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
%             xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
%             set (gca, 'fontsize', 20)
%             hold on 
%             plot(Points(1:3:end,1),Points(1:3:end,2),'r.', 'markersize',12)
%             hold on 
%             quiver(Points(1:3:end,1),Points(1:3:end,2), V0Z(1,1:3:end)', V0R(1,1:3:end)', 'b')
%             legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(1)), ' eV'), 'Location','northwest','Interpreter','latex');
%             axis equal
%         subplot(2,2,2)
%             hold on 
%             plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
%             ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
%             xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
%             set (gca, 'fontsize', 20)
%             hold on 
%             plot(Points(1:3:end,1),Points(1:3:end,2),'r.', 'markersize',12)
%             hold on 
%             quiver(Points(1:3:end,1),Points(1:3:end,2), V0Z(5,1:3:end)', V0R(5,1:3:end)', 'b')
%             legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(5)), ' eV'), 'Location','northwest','Interpreter','latex');
%             axis equal
%         subplot(2,2,3)
%             hold on 
%             plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
%             ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
%             xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
%             set (gca, 'fontsize', 20)
%             hold on 
%             plot(Points(1:3:end,1),Points(1:3:end,2),'r.', 'markersize',12)
%             hold on 
%             quiver(Points(1:3:end,1),Points(1:3:end,2), V0Z(10,1:3:end)', V0R(10,1:3:end)', 'b')
%             legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(10)), ' eV'), 'Location','northwest','Interpreter','latex');
%             axis equal
%         subplot(2,2,4)
%             hold on 
%             plot(cathode.coefs(:,1),cathode.coefs(:,2), 'k-', 'linewidth',2) % plot cathode geometry
%             ylabel('$R$ [m]', 'interpreter', 'latex','Fontsize', 22)
%             xlabel('$Z$ [m]', 'interpreter', 'latex', 'Fontsize', 22)
%             set (gca, 'fontsize', 20)
%             hold on 
%             plot(Points(1:3:end,1),Points(1:3:end,2),'r.', 'markersize',12)
%             quiver(Points(1:3:end,1),Points(1:3:end,2), V0Z(end,1:3:end)', V0R(end,1:3:end)', 'b')
%             legend('cathode', '$e^-$',strcat('$\mathbf{v_0}$',' : E =',num2str(E(end)), ' eV'), 'Location','northwest','Interpreter','latex');
%             axis equal
%             
        
