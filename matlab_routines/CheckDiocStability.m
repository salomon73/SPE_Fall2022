t2dlength=size(M.t2d,1);
fieldstart=max(1,t2dlength-500);%floor(0.95*t2dlength);
deltat=t2dlength-fieldstart;
maxl=50;
nbsave=200;

%[~, rgridend]=min(abs(M.rgrid-0.045));
rgridend=sum(M.nnr(1:2));
[~, name, ~] = fileparts(M.file);

dens=mean(M.N(:,:,fieldstart:end),3);
[R,Z]=meshgrid(M.rgrid,M.zgrid);
Rinv=1./R;
Rinv(isnan(Rinv))=0;
VTHET=mean(M.fluidUTHET(:,:,fieldstart:end),3);
omegare=(VTHET.*Rinv');
Er=mean(M.Er(:,:,fieldstart:end),3);
Ez=mean(M.Ez(:,:,fieldstart:end),3);
vdrift=(M.Br'.*Ez-Er.*M.Bz')./M.B'.^2;
omegadrift=(vdrift.*Rinv');
nbzid=9;
zindices=linspace(length(M.zgrid)/6,length(M.zgrid)*5/6,nbzid);

for zid=1:nbzid%1:length(zindices)%;%floor(length(M.zgrid)/2);
    zindex=floor(zindices(zid));
    rindex=M.geomweight(:,zindex,1)>=0;
    
    n=dens(rindex,zindex);
    omegar=omegadrift(rindex,zindex);
    rgrid=M.rgrid(rindex);
    Bz=mean(M.Bz(zindex,rindex),2);
    zfolder=sprintf('diocz_%03d',zindex);
    if ~isfolder(zfolder)
        mkdir(zfolder)
    end
    wce=abs(M.qsim/M.msim*Bz);
    wpe2=M.qsim^2/(M.msim*M.eps_0*M.weight)*n;
    wnum=complex(zeros(nbsave,maxl));
    rnumlength=max(length(rgrid)+2,4003);
    phinum=complex(zeros(nbsave,rnumlength,maxl));
    rnum=zeros(rnumlength,maxl);
    parfor l=1:maxl
        [w, phin, r] = diocotronstab(l, rgrid, n, omegar, wce, 1024+512);
        
        % different profiles of dphi
        [~,ind]=sort(abs(imag(w)),'descend');
        
        wnum(:,l)=w(ind(1:nbsave));
        phinum(:,:,l)=padarray(phin(:,ind(1:nbsave))',[0 rnumlength-size(phin,1)], 0, 'post');
        rnum(:,l)=padarray(r,[0 rnumlength-length(r(:))],0, 'post')';
        sprintf('l=%02d, zid=%02d done',l,zid)
    end
    lsteps=1:maxl;
    savecase(zfolder,rnum,wnum,phinum,rgrid,wce,Er,Ez,n,omegare,vdrift,omegadrift,M.geomweight,rgrid(1),rgrid(end),lsteps);
end

function savecase(folder,r,w,phinum,rgrid,wce,Er,Ez,n,omegare,vdrift,omegadrift,geomweight,a,b,lsteps)
    save(sprintf('%s/results',folder),'r','w','phinum','a','b','lsteps','rgrid','wce','Er','Ez','n','omegare', 'vdrift','omegadrift','geomweight','phinum')
end