function [Y] = care_nwt_fac(Y0,A,B,C,tol,maxsteps)
% Newton's method for continuous-time algebraic Riccati equations
%  (CARE)  0  =  C'C  +  A' X  +  X A  -  X BB' X  =: R(X)
%
% author  Peter Benner
% date    2005/12/23
% 
%
% Input:
%  Y0        initial starting guess s.t. A - BB'Y_0'Y_0 is stable.
%                  Note: this is not checked - if Y_0 is not stabilizing, then
%                  the iteration may fail to converge or converge to a
%                  non-stabilizing solution!
%  A         Matrix from (CARE)
%  B         Matrix from (CARE)
%  C         Matrix from (CARE)
%  tol       stopping criterion, i.e. the iteration ist stopped if
%                  |R(X)|_F/max(1,|Y'Y|_F) <= tol 
%                  , default sqrt(eps*n)
%
%  maxsteps  maximum number of iteration steps, default 50
%
% Output:
%  Y        approximate factor solution of CARE so that X=Y'Y

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
% Copyright (C) Peter Benner 2005
%

narginchk(4,6)
n = size(A,1);
m = size(B,2);
p = size(C,1);

%% check matrix sizes
if size(A,2) ~= n,  
    error('A must be square.'),  
end
if size(B,1) ~= n,  
    error('B must have the same number of rows as A.'),  
end
if size(C,2) ~= n,  
    error('C must have the same number of columns as A.'),  
end

if nargin < 6,  maxsteps = 50;  end
if (nargin < 5),  
  tol = sqrt(eps*n);
else
    if tol < sqrt(eps*n)
        tol = sqrt(eps*n);
        warning('Error tolerance too small, may not be achieved!')
    end
end

%% Initialization
iter = 0;
if isempty(Y0),
    Y = zeros(1,n);
else
    if size(Y0,2) ~= n,  
        error('Y0 must have the same number of rows as A.'),  
    else
        Y = Y0;
    end
end
YA    = Y*A;
YB    = Y*B;
nres  = norm(C'*C + YA'*Y + Y'*YA - YB*YB','fro');     
Xnorm = norm(Y*Y','fro');
Err = nres/max(1,Xnorm);
onemore     = 0;
convergence = Err <= tol;

%% Newton iteration
while (iter < maxsteps) && ((~convergence) || (convergence && (onemore < 2))),
  % Here one may employ RRQR to compress W.
  W        = [C; YB'*Y]; 
  Y        = lyap_sgn_fac(A - B*(YB)'*Y,W); 
  YA       = Y*A;
  YB       = Y*B;
  nres     = norm(C'*C + YA'*Y + Y'*YA - (Y'*YB)*(YB'*Y),'fro');     
  Xnorm    = norm(Y*Y','fro');
  iter = iter + 1;
% Uncomment next line for verbose mode.
 fprintf('||R(X_%i)||/||X|| = %d', iter, nres/Xnorm)
  Err = nres/max(1,Xnorm);
  convergence = Err <= tol;
  if convergence,  onemore = onemore + 1;  end
end

if (iter == maxsteps) && (nres/max(1,Xnorm) > tol),
  fprintf('CARE_NWT_FAC: no convergence in %d iterations', maxsteps)
end
