function [pot] = PotentialWell(Phi,Atheta,r,z,creategraph)
%PotentialWell Computes and display the potential well given the electrostatic potential Phi 
% and the magnetic vector potential Atheta
if nargin <5
    creategraph=true;
end
[rgrid,~]=meshgrid(r,z);
rAthet=(rgrid.*Atheta)';
contpoints=contourc(z,r,rAthet,rAthet(:,1)');
[x,y,zcont]=C2xyz(contpoints);
k=1;
zdiff=[diff(zcont),0];
for j=1:size(x,2)
    %for i=1:size(pot,1)
        xloc=x{j};
        yloc=y{j};
        xloc=xloc(yloc<r(end) & yloc>r(1));
        yloc=yloc(yloc<r(end) & yloc>r(1));
        x{j}=xloc;
        y{j}=yloc;
        if(length(xloc)>1)
            pot{j}=interp2(z,r,Phi,xloc,yloc);
            pot{j}=pot{j}-min(pot{j});
        end
        k=k+(zdiff(j)~=0);
end
x=cell2mat(x);
y=cell2mat(y);
pot=cell2mat(pot);
[Z,R]=meshgrid(z,r);
pot=griddata(x,y,pot,Z,R);
if(creategraph)
figure
surface(z(2:end),r(2:end),-pot(2:end,2:end),'edgecolor','none')
xlabel('r')
ylabel('z')
colorbar
end
end
