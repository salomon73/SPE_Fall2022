function [Continuity]=Continuityspline(M,it,norm,log,contourscale, zlims)
%ForceBalance Show the radial force balance
%   Plot the three radial forces for the given time-step it or averaged
%   over the range of time steps defined in it

if strcmp(it,':')
    it=floor(0.95*size(M.t2d)):size(M.t2d);
end
switch nargin
    case 2
        contourscale=0.2;
        norm=false;
        log=false;
        zlims=[-inf inf];
    case 3
        log=false;
        contourscale=0.2;
        zlims=[-inf inf];
    case 4
        contourscale=0.2;
        zlims=[-inf inf];
    case 5
        zlims=[-inf inf];
    case 6
        
    otherwise
        error("Invalid number of arguments")
end

contourcolor=[0 0 0];%[255 20 147]/255;

if(it(1)==1)
    it=it(2:end);
end
if(it(end)==length(M.t2d))
    it=it(1:end-1);
end

N=M.N(:,:,it);
maxdens=max(N(:));
densitycontour=contourc(M.zgrid,M.rgrid,mean(N,3),[contourscale contourscale]*maxdens);
rgridmax=sum(M.nnr(1:2));

%% Continuity
t=M.t2d(it);
dndt=zeros(size(N,1),size(N,2),size(N,3));
S=zeros(size(N,1),size(N,2),size(N,3));
ndivu=zeros(size(N,1),size(N,2),size(N,3),3);
ugradn=zeros(size(N,1),size(N,2),size(N,3),2);
[R,~]=meshgrid(M.rgrid,M.zgrid);
R=R';
Rinv=1./R;
Rinv(isinf(Rinv))=0;
geommask=M.geomweight(:,:,1)<0;
for i=1:size(N,3)
    ddt=(M.N(:,:,it(i)+1)-M.N(:,:,it(i)-1))./(M.t2d(it(i)+1)-M.t2d(it(i)-1));
    ddt(geommask)=0;
    
    dndt(:,:,i)=ddt;
    
    ndivu1=N(:,:,i).*(M.fluidUR.der(:,:,it(i),[1 0]));
    ndivu1(geommask)=0;
    ndivu(:,:,i,1)=ndivu1;
    
    ndivu2=N(:,:,i).*(M.fluidUR(:,:,it(i)).*Rinv);
    ndivu2(geommask)=0;
    ndivu(:,:,i,2)=ndivu2;
    
    ndivu3=N(:,:,i).*(M.fluidUZ.der(:,:,it(i),[0 1]));
    ndivu3(geommask)=0;
    ndivu(:,:,i,3)=ndivu3;
    
    ugradn1=  M.fluidUR(:,:,it(i)).*M.N.der(:,:,it(i),[1 0]);
    ugradn1(geommask)=0;
    ugradn(:,:,i,1)=ugradn1;
    ugradn2=  M.fluidUZ(:,:,it(i)).*M.N.der(:,:,it(i),[0 1]);
    ugradn2(geommask)=0;
    ugradn(:,:,i,2)=ugradn2;
    
    if M.neutcol.present && ~isempty(M.neutcol.io_cross_sec)
        % average kinetic energy in each direction
%         Er=M.fluidEkin(1,:,:,it(i));
%         Ethet=M.fluidEkin(2,:,:,it(i));
%         Ez=M.fluidEkin(3,:,:,it(i));
%         Ek=squeeze((Er+Ethet+Ez)/M.qe);
% %         % corresponding particle velocity
%         U=squeeze(sqrt(2*M.weight/M.msim*(Er+Ethet+Ez)));
         U=-M.Er(:,:,it(i))./M.Bz';
         Ek=0.5*M.msim/M.weight/M.qe*U.^2;
         Ek(N(:,:,i)<=0)=0;
        % ionisation cross section per cell
        sigio=M.sigio(Ek);
        Si=M.neutcol.neutdens.*N(:,:,i).*sigio.*U;
    end
    
    if M.maxwellsrce.present
        Si=Si+M.maxwellsrce.frequency*M.weight/(pi*diff(M.maxwellsrce.zlim)*(M.maxwellsrce.rlim(2)^2-M.maxwellsrce.rlim(1)^2));
    end
    Si(geommask)=0;
    S(:,:,i)=Si;
end
Continuity.N=N;
Continuity.dndt=dndt;
Continuity.ndivu=ndivu;
Continuity.ugradn=ugradn;
Continuity.maxdens=maxdens;
Continuity.it=it;

Continuity.S=S;
end

