function y=lyapunov(Z,x,eqn,oper,opts)
% Computes matrix vector product with the Lyapunov operator.
%
% Input:           
%  Z         Low-rank solution factor of the Riccati equation
%  x         vector for matrix vector product
%  eqn       structure with data for A, E and fields B, C, G  
%                  eqn.E(optional, eqn.haveE specifies whether it is
%                  there) in the above equation with ZZ' approximating X
%                  eqn.haveUV specfies whether is feedback there
%                  
%  oper      structure contains function handles for operations with
%                  A, E
%  opts      structure contains the following fields
%  eqn.type 'N' or 'T' for type of lyapunov equation
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

if eqn.type=='N'
  adjoint='T';
else
  adjoint='N';
end

if eqn.haveE
  z=Z*(Z'*(oper.mul_E(eqn, opts,adjoint,x,'N')));
  if eqn.haveUV
    y = oper.mul_A(eqn, opts,eqn.type,z,'N') - ...
      eqn.V*(eqn.U'*z)           + ...
      oper.mul_E(eqn, opts,eqn.type, Z*Z'*(oper.mul_A(eqn, opts,adjoint,x,'N')) , 'N');
  else             % No feedback
    y= oper.mul_A(eqn, opts,eqn.type,z,'N') + ...
      oper.mul_E(eqn, opts,eqn.type,Z*(Z'*(oper.mul_A(eqn, opts,adjoint,x,'N'))), 'N');
  end
  if isfield(opts,'nm')
    y=y+eqn.G*(eqn.G'*x);
  else
    if eqn.type=='N'
      y=y+eqn.B*(eqn.B'*x);
    else
      y=y+eqn.C'*(eqn.C*x);
    end
  end
else               % Standard equation
  z=Z*(Z'*x);
  if eqn.haveUV    
    y=  oper.mul_A(eqn, opts,eqn.type,z,'N') - ...
      eqn.V*(eqn.U'*z)+ Z*(Z'*(oper.mul_A(eqn, opts,adjoint,x,'N')));
  else             % No feedback
    y=  oper.mul_A(eqn, opts,'N',z,'N') + ...
      Z*(Z'*(oper.mul_A(eqn, opts,adjoint,x,'N')));
  end
  if isfield(opts,'nm')
    y=y+eqn.G*(eqn.G'*x);
  else
    if eqn.type=='N'
        y=y+eqn.B*(eqn.B'*x);
    else
      y=y+eqn.C'*(eqn.C*x);
    end
  end
end

% in case of Rosenbrock we get a -1/(2*timestep)*(E'*Z*Z'*E)
% from both F and F'
if isfield(opts,'rosenbrock')&&~isempty(opts.rosenbrock)
  if eqn.haveE       % generalized equations
    y=y-(1/opts.rosenbrock.stepsize)*oper.mul_E(eqn, opts,eqn.type,z,'N');
  else
    y=y-(1/opts.rosenbrock.stepsize)*z;
  end
end
