% ~~~~~~~~~~~~~~~~~~~~~~~~
function dXdt = twobodyEOMJ2(t,X)
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

r = norm(rvec) ;
rz = rvec(3);
% J2 terms
J2Terms =  (1.5*J2*(Re/r)^2*...
    [1-5*rz^2/r^2;...
     1-5*rz^2/r^2;...
     3-5*rz^2/r^2]);

rdotvec = vvec ;
vdotvec = -(mu/r^3)*rvec .* (1 + J2Terms);

dXdt = [rdotvec; vdotvec];
end 