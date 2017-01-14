function X =mess_dense_nm(A,B,C,E,X0)
% naive Newton Kleinman iteration  for the ARE 
%  
%     A'*X*E+E'*X*A+C'*C-E'*X*B*B'*X*E = 0
%
% Inputs:
%
% A,B,C,E    Coefficients in the above equation
% X0         initial guess for the solution
%
% Outputs:
% X          Solution 
%  

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
  
% Author: Jens Saak
  tol = 1e-12;
  maxiter = 50;
  
  F = B*B';
  G = C'*C;
  
  res0 = norm(G);
  
  if (nargin<4) || isempty(E)
    E = eye(size(A,1));
  end
  for i=1:maxiter
    if i>1 || ( (nargin==5) && ~isempty(X0) )
      K = B'*X0*E;
    else
      K = zeros(size(B'));
      X0 = zeros(size(A));
    end
    X = lyap(A-K'*B',G+K'*K,[],E);
    XE = X*E;
    res = norm(A'*XE+XE'*A-XE'*F*XE+G);
    rc = norm(X-X0)/norm(X);
    if (rc<tol) || (res<tol*res0)
      break
    else
      X0 = X;
    end
  end
