% ~~~~~~~~~~~~~~~~~~~~~~~~
function dXdt = twobodyEOMJ2WithSTM(t,X)
% ------------------------
%{
  This function calculates first and second time derivatives of r
  governed by the equation of two-body 2D motion with j2 pertubations for
  satellites orbiting Earth
 
  X    - column vector containing r and v at time t
  dXdt - column vector containing drdt and dvdt at time t

  User M-functions required: none
%}
% ~~~~~~~~~~~~~~~~~~~~~~~~

mu = getConstants().Earth.mu;
J2 = getConstants().Earth.J2;
Re = getConstants().Earth.Re;
rvec = X(1:3);
vvec = X(4:6);
Phivec = X(7:end);
Phi = [reshape(Phivec,6,6)];

r = norm(rvec) ;
rx = rvec(1);
ry = rvec(2);
rz = rvec(3);
% J2 terms
J2Terms =  (1.5*J2*(Re/r)^2*...
    [1-5*rz^2/r^2;...
     1-5*rz^2/r^2;...
     3-5*rz^2/r^2]);

rdotvec = vvec ;
vdotvec = -(mu/r^3)*rvec .* (1 + J2Terms);

dv_dr = [(-(mu/r^3)+3*mu*rx^2/r^5)*(1 + J2Terms(1))-(mu*rx/r^3)*((-3*J2*Re^2*rx/r^4)*(1-5*rz^2/r^2)+(1.5*J2*(Re/r)^2)*(10*rx*rz^2/r^4)),...
         (3*mu*rx*ry/r^5)*(1 + J2Terms(1))-(mu*rx/r^3)*((-3*J2*Re^2*ry/r^4)*(1-5*rz^2/r^2)+(1.5*J2*(Re/r)^2)*(10*ry*rz^2/r^4)),...
         (3*mu*rx*rz/r^5)*(1 + J2Terms(1))-(mu*rx/r^3)*((-3*J2*Re^2*rz/r^4)*(1-5*rz^2/r^2)-(1.5*J2*(Re/r)^2)*(10*rz*(ry^2+rx^2)/r^4));
         ...
         (3*mu*rx*ry/r^5)*(1 + J2Terms(2))-(mu*ry/r^3)*((-3*J2*Re^2*rx/r^4)*(1-5*rz^2/r^2)+(1.5*J2*(Re/r)^2)*(10*rx*rz^2/r^4)),...
         (-(mu/r^3)+3*mu*ry^2/r^5)*(1 + J2Terms(2))-(mu*ry/r^3)*((-3*J2*Re^2*ry/r^4)*(1-5*rz^2/r^2)+(1.5*J2*(Re/r)^2)*(10*ry*rz^2/r^4)),...
         (3*mu*ry*rz/r^5)*(1 + J2Terms(2))-(mu*ry/r^3)*((-3*J2*Re^2*rz/r^4)*(1-5*rz^2/r^2)-(1.5*J2*(Re/r)^2)*(10*rz*(ry^2+rx^2)/r^4));...
         ...
         (3*mu*rx*rz/r^5)*(1 + J2Terms(3))-(mu*rz/r^3)*((-3*J2*Re^2*rx/r^4)*(3-5*rz^2/r^2)+(1.5*J2*(Re/r)^2)*(10*rx*rz^2/r^4)),...
         (3*mu*ry*rz/r^5)*(1 + J2Terms(3))-(mu*rz/r^3)*((-3*J2*Re^2*ry/r^4)*(3-5*rz^2/r^2)+(1.5*J2*(Re/r)^2)*(10*ry*rz^2/r^4)),...
         (-(mu/r^3)+3*mu*rz^2/r^5)*(1 + J2Terms(3))-(mu*rz/r^3)*((-3*J2*Re^2*rz/r^4)*(3-5*rz^2/r^2)-(1.5*J2*(Re/r)^2)*(10*rz*(ry^2+rx^2)/r^4))];

A = [zeros(3), eye(3);...
    dv_dr, zeros(3)];
Phidot = A*Phi;
Phidotvec = reshape(Phidot, 6*6,1);
dXdt = [rdotvec; vdotvec; Phidotvec];
end 