function [Phi,Er] = Analytic_stpe_dens_Efield(a,b,deltaphi,rm,rp,n,r)
%Analytic_stpe_dens_Efield Computes the electric field and potential for a
% step density in a biased coaxial geometry
%

q=1.602176620000000e-19;
eps0=8.854187817620000e-12;
phia=-deltaphi;
phib=0;


phip=1/log(b/a)*(q*n/2/eps0*rm^2*log(b/rp)*log(rp/rm)+(phia + q*n/2/eps0*(rp^2-rm^2)*(log(rp/a)-0.5))*log(b/rp));
phim=1/log(b/a)*(-q*n/2/eps0*rm^2*log(rm/a)*log(rp/rm)+( q*n/2/eps0*(rp^2-rm^2)*(log(b/rp)+0.5))*log(rm/a)+phia*log(b/rm));

Phi=zeros(size(r));
Er=zeros(size(r));

Phi(r>=a & r<rm)=((phim-phia)*log(r(r>=a & r<rm))+phia*log(rm)-phim*log(a))/log(rm/a);
Phi(r>=rm & r<rp)=-q*n/4/eps0*(r(r>=rm & r<rp).^2-rp^2) ...
    +(phip-phim +q*n/4/eps0*(rp^2-rm^2))/log(rp/rm)*log(r(r>=rm & r<rp)/rp)...
    +phip;
Phi(r>=rp & r<=b)=((phib-phip)*log(r(r>=rp & r<=b))+phip*log(b)-phib*log(rp))/log(b/rp);

Er(r>=a & r<rm)=-(phim-phia)./r(r>=a & r<rm)/log(rm/a);
Er(r>=rm & r<rp)=+q*n/2/eps0*r(r>=rm & r<rp) ...
    -(phip-phim +q*n/4/eps0*(rp^2-rm^2))/log(rp/rm)./r(r>=rm & r<rp);
Er(r>=rp & r<=b)=-(phib-phip)./r(r>=rp & r<=b)/log(b/rp);

end

