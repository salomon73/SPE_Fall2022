function out = PartInfos(nElectrons, nEnergy, nComponents)

    PartInfoV0 = zeros(5,npartsV0);
    PartInfoVn = zeros(5,nparts);
    compteur   = 0;
    for ii = 1:length(E)

       for jj = 1: nElectrons

           for kk = 1: nComponents
               
                compteur = compteur + 1; 
                PartInfoV0(1,compteur) = compteur;  
                PartInfoV0(2,compteur) = Points(1,jj);
                PartInfoV0(3,compteur) = Points(2,jj);
                PartInfoV0(4,compteur) = V0(ii,jj,kk,1);
                PartInfoV0(5,compteur) = V0(ii,jj,kk,2);
           end

       end

    end

    compteur = 0;
    for ii =1:length(E)

        for jj =1: nElectrons
            
                compteur = compteur + 1; 
                PartInfoVn(1,compteur) = compteur;  
                PartInfoVn(2,compteur) = Points(1,jj);
                PartInfoVn(3,compteur) = Points(2,jj);
                PartInfoVn(4,compteur) = VNorm(ii,jj,1);
                PartInfoVn(5,compteur) = VNorm(ii,jj,2);
            

        end
    end
end