part=M.species(1);
time=part.tpart;
track=true;
gcs=true;
rpos(1,1,:)=M.rgrid;
zpos(1,1,:)=M.zgrid;
partindices=M.species(1).partindex(1:100,end);
partindices=partindices(partindices>0);
timeindices=(-200:1:0)+length(part.tpart);
timeindices(timeindices<1)=[];
pR=part.R(partindices,timeindices,track);
pZ=part.Z(partindices,timeindices,track);
[~,Rind]=min(abs(pR-rpos),[],3);
[~,Zind]=min(abs(pZ-zpos),[],3);
locB=interp2(M.zgrid,M.rgrid,M.B',pZ,pR,'makima');
Vpar=part.Vpar(partindices,timeindices,track,gcs);
Vperpst=part.Vperp(partindices,timeindices,track,gcs);
must=0.5*M.me*Vperpst.^2./locB;
mustm=zeros(size(must));
for i=1:size(must,1)
mustm(i,:)=must(i,:)/mean(must(i,:));
end

lnstyleorder={'-','--','-.',':'};
figure();
semilogy(time(timeindices),mustm)
ax=gca;
ax.LineStyleOrder=lnstyleorder;
xlabel('time')
ylabel('mu*')

Vperp=part.Vperp(partindices,timeindices,track,false);
mu=0.5*M.me*Vperp.^2./locB;
lnstyleorder={'-','--','-.',':'};
mum=zeros(size(mu));
for i=1:size(mu,1)
mum(i,:)=mu(i,:)/mean(mu(i,:));
end
figure();
semilogy(time(timeindices),mum)
ax=gca;
ax.LineStyleOrder=lnstyleorder;
xlabel('time')
ylabel('mu')

Phi=zeros(size(Rind));
for i=1:size(Rind,2)
    [~,tfield]=min(abs(M.t2d-time(i)));
    timPhi=M.pot(:,:,tfield);
    %posindPhi=sub2ind(size(timPhi),Rind(:,i),Zind(:,i));
    Phi(:,i)=interp2(M.zgrid,M.rgrid,timPhi,pZ(:,i),pR(:,i));%timPhi(posindPhi);
end
Ebar=0.5*M.me*Vpar.^2+0.5*M.me*Vperpst.^2-M.qe*Phi;
ebarm=zeros(size(Ebar));
for i=1:size(ebarm,1)
ebarm(i,:)=Ebar(i,:)/mean(Ebar(i,:));
end
figure();
semilogy(time(timeindices),ebarm)
ax=gca;
ax.LineStyleOrder=lnstyleorder;
xlabel('time')
ylabel('Ebar/<Ebar>_t')

Phi=zeros(size(Rind));
for i=1:size(Rind,2)
    [~,tfield]=min(abs(M.t2d-time(i)));
    timPhi=M.pot(:,:,tfield);
    %posindPhi=sub2ind(size(timPhi),Rind(:,i),Zind(:,i));
    Phi(:,i)=interp2(M.zgrid,M.rgrid,timPhi,pZ(:,i),pR(:,i));%Phi(:,i)=timPhi(posindPhi);
end
E=0.5*M.me*Vpar.^2+0.5*M.me*Vperp.^2-M.qe*Phi;
em=zeros(size(E));
for i=1:size(em,1)
em(i,:)=E(i,:)/mean(E(i,:));
end
figure();
semilogy(time(timeindices),em)
ax=gca;
ax.LineStyleOrder=lnstyleorder;
xlabel('time')
ylabel('E/<E>_t')

R=pR';
Z=pZ';
Thet=part.THET(partindices,timeindices,true)';
figure();
plot(Z,R)
ax=gca;
ax.LineStyleOrder=lnstyleorder;
xlabel('Z')
ylabel('R')

figure();
plot(R.*cos(Thet),R.*sin(Thet))
ax=gca;
ax.LineStyleOrder=lnstyleorder;
xlabel('X')
ylabel('Y')