function [w, phin, r] = diocotronstab(l, R, n, wre, wce, nRmin)
if nargin<6
    nRmin=512;
end
% Necessary physical constants
% electron charge and mass
    e = 1.602176634e-19;
    me = 9.1093837015e-31;
% vacuum permittivity
    eps_0 = 8.854188e-12; 

% declaration of wpe2 (vector for numerical solution, renormalised on wce)
   
    wpe2 = n*(e^2)/(me*eps_0);
    
% scale factor for the frequency
    wscale=max(abs(wpe2/(2*wce)));
    
 % construction of r, new vector of distances (chooses the thinner mesh and
 % interpolate all variables on that mesh
   
    delta = diff(R);
    delta=delta(delta>0);
    x = unique(delta);
    N = numel(x);
    count = zeros(N,1);
    for k = 1:N
        count(k) = sum(delta==x(k));
    end
    x(count<2)=[];
    rmin = min(x);
    nR = floor((max(R)-min(R))/rmin) + 1;
    if nR<nRmin
        nR=nRmin+1;
    end
    r = linspace(R(1),R(end),nR);
    b=max(R);
    rmin=rmin/b;
    r=r/b;
    R=R/b;
 
 % interpolation  
 % of n, wre and wpe2 --> nin, wrin, wpe2in
   wrin = interp1(R(2:end),wre(2:end),r,'spline');
   wpe2in = interp1(R,wpe2,r);
   wpe2in=wpe2in/wscale^2;
   wrin=wrin/wscale;
   wce=wce/wscale;
   
 % finite difference to compute d/dr(wpe^2)
    Dwpe2 = zeros([nR 1]);
    Dwpe2(2:(nR-1)) = (wpe2in(3:end)-wpe2in(1:(nR-2)))/(2*rmin);
    
%  construction of matrices A & B such that A*phi = w*B*phi
    
    A = zeros(nR,nR);
    B = zeros(nR,nR);

    
   %    finite differences in A & B + other terms of the eigenvalue
   %    equation
    
    for i = 2:1:(nR-1)
            A(i,i) =   l*wrin(i)*(-2/(rmin)^2 - l^2/r(i)^2) - l*Dwpe2(i)/wce/r(i);
            A(i,i-1) = l*wrin(i)*(1/(rmin)^2 - 1/(2*r(i)*rmin));
            A(i,i+1) = l*wrin(i)*(1/(rmin)^2 + 1/(2*r(i)*rmin));
            
            B(i,i) = - 2/(rmin)^2 - (l^2/r(i)^2);
            B(i,i-1) = 1/(rmin)^2 - 1/(2*r(i)*rmin);
            B(i,i+1) = 1/(rmin)^2 + 1/(2*r(i)*rmin);
    
    end
    A=A(2:(end-1),2:(end-1));
    B=B(2:(end-1),2:(end-1));
    
% Solution of the eigenvalue problem
    [phin,w] = eig(A, B); 
 
% renormalisation of w and r
    w=diag(w);
    w = w*wscale;
    r=r*b;
end

