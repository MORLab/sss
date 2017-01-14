function y=riccati(Z,x,eqn,oper,opts)
% Computes matrix vector product with the Riccati operator.
%  y = A^T*ZZ^T*x + ZZ^T*Ax + C^T*Cx - ZZ^T*BB^T*ZZ^T*x  or    
%  y = A^T*ZZ^T*Ex + C^T*Cx + E^T*ZZ^T*Ax - BB^T*ZZ^Tx  
%  
% Input:                
%  Z         Low-rank solution factor of the Riccati equation
%  x         vector for matrix vector product
%  eqn       structure with data for A, E and fields B, C  
%                  eqn.E(optional, eqn.haveE specifies whether it is
%                  there) in the above equation with ZZ' approximating X
%  oper      structure contains function handles for operations with
%                  A, E
%  opts    structure contains parameters for the algorithm
%
% Output:
%  y        result of matrix vector product

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

% generalized equations
% y = A'*Z*Z'*E*x + C'*C*x + E'*(Z*Z'*A*x - B*B'*Z*Z'*x)  
%
% uses operatorfunctions mul_E, mul_A

%% Check input
if ~isfield(opts,'bdf'), opts.bdf=[]; end
if isstruct(opts.bdf) && isfield(opts.bdf, 'tau') && isfield(opts.bdf, 'beta')
    bdf = 1;
    pc = -1 / (2 * opts.bdf.tau * opts.bdf.beta);
else
    bdf = 0;
end

%% Compute MVP
if eqn.haveE       
  z=Z*(Z'*(oper.mul_E(eqn, opts,'N',x,'N')));
  if bdf
      y1 = (opts.bdf.tau * opts.bdf.beta) * ...
          oper.mul_ApE(eqn, opts,'T',pc, 'T', z,'N') + eqn.C'*(eqn.C*x);
      y2 = Z*(Z'*((opts.bdf.tau * opts.bdf.beta) * ...
          oper.mul_ApE(eqn, opts,'N',pc, 'N', x,'N')-eqn.B*(eqn.B'*z)));
  elseif eqn.setUV                                                        
      y1 = oper.mul_A( eqn, opts, 'T', z, 'N' ) ...                   
          - eqn.V(:,1:eqn.UVsize) * (eqn.U(:,1:eqn.UVsize)' * z) ... 
          + eqn.C' * (eqn.C * x);                                    
      y2 = Z * (Z' * (oper.mul_A( eqn, opts, 'N', x, 'N' ) ...        
          - eqn.U(:,1:eqn.UVsize) * (eqn.V(:,1:eqn.UVsize)' * x) ... 
          - eqn.B * (eqn.B' * z)));
  else
      y1 = oper.mul_A(eqn, opts,'T',z,'N') + eqn.C'*(eqn.C*x);
      y2 = Z*(Z'*(oper.mul_A(eqn, opts,'N',x,'N')-eqn.B*(eqn.B'*z)));
  end
  y  = y1 + oper.mul_E(eqn, opts,'T',y2,'N');
   
else
  z=Z*(Z'*x);
  if bdf
      y= (opts.bdf.tau * opts.bdf.beta) * ...
          oper.mul_ApE(eqn, opts,'T',pc, 'T', z,'N')        + ...
          Z*(Z'*((opts.bdf.tau * opts.bdf.beta) * ...
          oper.mul_ApE(eqn, opts,'N',pc, 'N',x,'N') - ...
          eqn.B*eqn.B'*z))                + ...
          eqn.C'*(eqn.C*x);
  elseif eqn.setUV                                                        
      y = oper.mul_A( eqn, opts, 'T', z, 'N' ) ...                    
          - eqn.V(:,1:eqn.UVsize) * (eqn.U(:,1:eqn.UVsize)' * z) ...  
          + Z * (Z' * (oper.mul_A( eqn, opts, 'N', x, 'N' ) ...       
          - eqn.U(:,1:eqn.UVsize) * (eqn.V(:,1:eqn.UVsize)' * x) ...  
          - eqn.B * eqn.B' * z)) ...                                  
          + eqn.C' * (eqn.C * x);
  else
      y= oper.mul_A(eqn, opts,'T',z,'N')        + ...
          Z*(Z'*(oper.mul_A(eqn, opts,'N',x,'N') - ...
          eqn.B*eqn.B'*z))                + ...
          eqn.C'*(eqn.C*x);
  end
end
