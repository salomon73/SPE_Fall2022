rbm=0.005;
rbp=0.0128;
%rbm=0.09;
%rbp=0.113;
%b=0.16;
%a=0.001;
a=M.rgrid(1);
b=M.rgrid(end);

Q=mean(M.Er(1,floor(length(M.zgrid)/2),end-30:end),3)/2*a;
%Q=mean(M.Er(1,100,end-30:end),3)/2*a;
omega_d=M.omepe^2/(2*M.omece);
if( a~=0)
    omega_q=2*Q/(M.B0*(rbm^2));
else
    omega_q=0;    
end

Levaluated=10;
stability=zeros(Levaluated,1);
for l=1:Levaluated
    stability(l)=(-l*(1-omega_q/omega_d)*(1-(rbm/rbp)^2)*(1-(a/b)^(2*l))+ ...
        2*(1+(a/b)^(2*l))-(1+(rbm/rbp)^(2*l))*((a/rbm)^(2*l)+(rbp/b)^(2*l)))^2 -...
        4*(rbm/rbp)^(2*l)*(1-(rbp/b)^(2*l))^2*(1-(a/rbm)^(2*l))^2;
    if(stability(l) <0)
        fprintf(2,"Diocotron mode l=%d is potentially unstable\n",l);
    end
end
fig=figure;
plot(1:Levaluated,stability)
hold on
plot([1 Levaluated],[0 0],'k')
for l=1:Levaluated
    if stability(l)<0
        plot(l,stability(l),'rx')
    end
end
xlabel('azimuthal mode number')
ylabel('stability condition')
title(sprintf('r_b^-=%1.2g [m] r_b^+=%1.2g [m] a=%1.2g [m] b=%1.2g [m]',rbm,rbp,a,b))
M.savegraph(fig,sprintf('%s/%s_diocstab',M.folder,M.name),[15,10])