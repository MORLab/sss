function [M,E,K]=genMEK(n)
%  [M,E,K]=genMEK(n)
%
%  Generate system matrices of a mass-spring-damper-system 
%
%              M x''(t)  = -E x'(t) -Kx(t)
%
%  for random spring, mass and damping data.
%
%  INPUT:
%  n   desired dimension of the resulting system matrices
%
%  OUTPUT:
%  M   mass matrix of the system, i.e., diagonal matrix containing
%      the masses.
%
%  E   damper matrix containing the damping coefficients in the
%      proportional damping.
%
%  K   stiffness matrix, i.e., tridiagonal matrix built from the
%      stiffnesses of the springs.

M=spdiags(rand(n,1),0,n,n);
while (any(diag(M)==0))
  M=spdiags(rand(n,1),0,n,n);
end

E=spdiags(1e-2*rand(n,1),0,n,n);
while (any(diag(E)<=0))
  E=spdiags(rand(n,1),0,n,n);
end

%x=rand(n,1);
x=ones(n,1);
y=[x(1:n-1)+x(2:n); x(n)];
z=[x(n); x(1:n-1)];

K= spdiags( [-x y -z], -1:1,n,n);