function [nrm,k,T,V] = res2_norms(Z,Rmul,eqn,fopts,oper)
% Computes the 2 Norm of the residual of Z for the symmetric operator given 
%   by y=Rmul(Z, x, eqn, opts,oper);
% e.g., for the generalized Lyapunov residual Rmul implements
%   R := FZZ^T*E^T + EZZ^T*F^T + GG^T (1) 
% 
%
% That means res2 computes the spectral radius of this operator R by a
% Lanzcos iteration. The function Rmul should exploite the structure of F 
% and rectangular structure of Z and G. Thus it can be computed in O(n) 
% effort and is therefore much cheaper than the computation of the e.g. the
% Frobenius norm.
%
% Method:
%
%    The Lanczos method produces matrices V and T such that
%
%      V(:,1) in span{rv},
%      V'*V = eye(k+1),
%      F*V(:,1:k) = V*T.
%
%  Remark:
%    V is only constructed if its required by the output arguments, or if
%    reorthogonalization is used.
% 
%    This implementation does not check for (near-)breakdown!
%
%  
% Input:                 
%  Z            Low-rank solution factor of operator
%
%  eqn          structure with data for operator 
%
%  Rmul         function handle to a function implementing the 
%                   multiplication with the residual operator of interest.                
%
%  opts         options structure with fields
%  opts.maxiter maximal number of Arnoldi steps (usually k<<n)
%                   (optional - chosen as 50 if omitted)
%  opts.tol     relative accuracy of the norm computation       
%                   (optional - chosen as 1e-6 if omitted)
%  opts.rv      initial n-vector
%                  (optional - chosen by random, if omitted)
%  opts.orth    reorthogonalization flag
%                  (optional - switched off, if omitted)
%
%  oper         structure contains function handles for operations with
%                  data in eqn
%
% Output:
%  nrm          the residual of the iterate X=Z*Z'
%
%  k            number of Lanzcos steps taken
%
%  T            matrix T ((k+1)-x-k matrix, symmetric tridiagonal);
%
%  V            matrix V (n-x-(k+1) matrix, orthogonal columns).
%
% uses eventually operatorfunctions in Rmul

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
%

%% check Rmul
n = size(Z,1); % Get system order.
switch Rmul
  case 'lyapunov'
    Rmul=@lyapunov;
    opts=fopts.adi;
  case 'riccati'
    Rmul=@riccati;
    opts=fopts.nm;
    otherwise
        errror('MESS:control_data','Rmul has to be lyapunov or riccati');
end

if ~isfield(opts,'res')
  error('MESS:control_data','residual control structure opts.res missing.');
end

if ~isfield(opts.res,'maxiter')||isempty(opts.res.maxiter)
    warning('MESS:control_data','opts.res.maxiter is set to 10 (default)');
    opts.res.maxiter = 10;
else
    if opts.res.maxiter >= n-1,
        error('maxiter must be smaller than the order of A!');
    end
end

if ~isfield(opts.res,'tol')||isempty(opts.res.tol),
    warning('MESS:control_data','opts.res.tol is set to 1e-6 (default)');
    opts.res.tol= 1e-6; 
end


if ~isfield(opts.res,'rv')||isempty(opts.res.rv),
    opts.res.rv = randn(n,1); 
end

if ~isfield(opts.res,'orth')||isempty(opts.res.orth),
    warning('MESS:control_data','opts.res.orth is set to 0 (default)');
    opts.res.orth=0; 
end

%% start computation
T = zeros(opts.res.maxiter+1,opts.res.maxiter);
v1 = (1.0/norm(opts.res.rv))*opts.res.rv;

% initial step
% Matrix-vector product R*v1
w=Rmul(Z,v1,eqn,oper,fopts);

T(1,1)=v1'*w;
r=w-T(1,1)*v1;
T(1,2)=norm(r);
v2=r./T(1,2);
T(2,1)=T(1,2)';

% store vectors in V if needed
if opts.res.orth || nargout==4
    V = zeros(n,opts.res.maxiter+1);
    V(:,1) = v1;
    V(:,2)=v2;
end

nrm = 0;

for k=2:opts.res.maxiter-1;
    %Lanczos 3-term recursion 
    
    % Matrix-vector product R*v2
    w=Rmul(Z,v2,eqn,oper,fopts);
    
    T(k,k)=v2'*w;
    r=w-T(k-1,k)*v1-T(k,k)*v2;
    
    %re-orthogonalization by MGS
    if opts.res.orth
        for j=1:k,
            r = r - (V(:,j)'*r)*V(:,j);
        end
    end
    
    T(k,k+1)=norm(r);
    v1=v2;
    v2=r./T(k,k+1);
    if opts.res.orth || nargout==4
        V(:,k+1)=v2;
    end
    T(k+1,k)=T(k,k+1)';
    nrmold= nrm;
    nrm = max(abs(eig(T(1:k,1:k))));
    if abs(nrm-nrmold)<opts.res.tol*nrm
        break;
    end
end
