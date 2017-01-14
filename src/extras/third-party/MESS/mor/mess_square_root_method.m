function [Er,Ar,Br,Cr,Dr,TL,TR] = mess_square_root_method(eqn,opts,oper,ZB,ZC)
% Square root method for the computation of the balanced and reduced
% system
%  
% Call
%  [Er,Ar,Br,Cr,Dr] = mess_square_root_method(eqn,opts,oper,ZB,ZC);
%
% Inputs:
%  eqn, opt, oper   the standard structures
%  ZB, ZC           the (tall and skinny) Gramian factors
%
% Outputs: 
%  Er,Ar,Br,Cr,Dr   The reduced system matrices.
%                   Note that Er is always an identity by
%                   construction and Dr will only be set for
%                   certain DAE systems.
%
% The implementation (especially in the case E!=I) follows the
% derivation in:
%	Efficient Numerical Solution of Large Scale Algebraic Matrix
%	Equations in PDE Control and Model Order Reduction; 
%   Saak, Jens;
%   Dissertation, TU Chemnitz; 2009. 
% 
% NOTE: Currently only standard state space systems and descriptor systems
% with E invertible are supported.
%

% Author: 
%  Jens Saak, Feb 24 2016

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

  [U0,S0,V0] = svd(ZC'*oper.mul_E(eqn, opts, 'N', ZB, 'N'),0);
  s0=diag(S0);
  ks=length(s0);
  k=ks;
  while (sum(s0(k-1:ks))<opts.bt.srm.tol/2)&&(k>2)
    k=k-1; 
  end
  k0=k;

  r= min([opts.bt.srm.max_ord k0]);
  if opts.bt.srm.info>0
    fprintf(1,'reduced system order: %d\n\n',r);
  end
  
  sigma_r=diag(S0(1:r,1:r));
  
  VB = ZB*V0(:,1:r);
  VC = ZC*U0(:,1:r);
  
  TR = VB*diag(ones(r,1)./sqrt(sigma_r));
  TL = VC*diag(ones(r,1)./sqrt(sigma_r));
  
  Ar = TL'*oper.mul_A(eqn, opts, 'N', TR, 'N');
  Br = TL'*eqn.B;
  Cr = eqn.C*TR;
  Er = eye(r);
  Dr = [];