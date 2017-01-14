function [Er,Ar,Br,Cr,S] = mess_tangential_irka(E,A,B,C,r,maxiter,tol,verbose)
% The tangential IRKA method with automatic selection of initial shifts and
% tangential directions.
%
% Call
%   [Er,Ar,Br,Cr,S] = mess_tangential_irka(E,A,B,C,r,maxiter,tol)
%
% Inputs:
%  E,A,B,C    The mass, system, input and output matrices describing the
%             original system
%  r          The reduced order
%  maxiter    maximum iteration number for the IRKA iteration
%  tol        bound for the relative change of the IRKA shifts used as
%             stopping criterion (usually 1e-2 to 1e-4 is sufficient)
%  verbose    1 : compute and show the sigma and error plots
%             0 : skip plotting (default)
%
% Outputs:
%  Er,Ar,Br,Cr,Dr   The reduced system matrices.
%
%
% NOTE: Currently only standard state space systems and descriptor systems
% with E invertible are supported.
%

% Author: 
%  Jens Saak, April 2016

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

if nargin<8, verbose=0; end

oper = operatormanager('default');
eqn.A_ = A;
eqn.E_ = E;
eqn.B = B;
eqn.C = C;
opts.dummy = 0;

[eqn, opts, oper] = oper.mul_A_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_pre(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);

%% Initialization
n = oper.size(eqn, opts);
m = size(B,2);
p = size(C,1);
% Choose the initial shifts as a set of random Ritz-values for (A,E)
U = orth(randn(n,r));
S = eig(full(U'*(oper.mul_A(eqn, opts, 'N', U, 'N'))),full(U'*(oper.mul_E(eqn, opts, 'N', U, 'N'))));
% and some random tangential directions
b = randn(m,r);
c = randn(p,r);

%% Start iteration
for iter = 1:maxiter
    S_old = S;
    %% Compute projection subspaces
    V = zeros(n,r);
    W = zeros(n,r);
    j = 1;
    while(j < r+1)        
        x = oper.sol_ApE(eqn, opts,'N',-S(j),'N',B*b(:,j),'N');
        y = oper.sol_ApE(eqn, opts,'T',-S(j),'T',C'*c(:,j),'N');
        if(abs(imag(S(j))) > 0)
            V(:,j) = real(x);
            W(:,j) = real(y);
            V(:,j+1) = imag(x);
            W(:,j+1) = imag(y);
            j = j + 2;
        else
            V(:,j) = real(x);
            W(:,j) = real(y);
            j = j + 1;
        end
    end    
    [V,~] = qr(V,0);
    [W,~] = qr(W,0);
   
    %% Compute ROM
    Er = (W'*oper.mul_E(eqn, opts, 'N', V, 'N'));
    Ar = (W'*oper.mul_A(eqn, opts, 'N', V, 'N'));
    Br = (W'*B);
    Cr = C*V;    
   
    %% Update interpolation points/tangential directions
    [T,S] = eig(Ar,Er);
    S = cplxpair(- diag(S), 1000*eps); % increased default pairing
                                       % tolerance since Octave
                                       % seems to give to
                                       % inaccurate egenvalues 
    b = (T\(Er\Br)).';
    c = Cr*T;
    
    %% Check for convergence
    err = norm(sort(S)-sort(S_old))/norm(S_old);
    %% comment out for non-verbose mode %%
    fprintf('IRKA step %3d, conv. crit. = %e \n', iter, err)

    if(err < tol)        
        break
    end
end
if ((iter == maxiter) && (err > tol)),
  warning('MESS:convergence','IRKA: No convergence in %d iterations.\n', maxiter)
end

if verbose
    mess_sigma_plot(eqn,opts,oper,Er,Ar,Br,Cr,[],-8,8,100);
end

[eqn, opts, oper] = oper.mul_A_post(eqn, opts, oper);
[eqn, opts, oper] = oper.mul_E_post(eqn, opts, oper);
[eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
