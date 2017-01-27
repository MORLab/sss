function [ V, eqn, opts, oper ]=mess_solve_shifted_system_BDF(eqn, opts, oper, pc, W)
% Solves (Ã + p*E)V = W for V, Ã = tau*beta*A - 0.5*E 
%  or Ã = tau*beta*A - 0.5*E - UV^T (BDF scheme)
%
% author  Björn Baran
% date    2015/09/01
%
%  Solves (Ã + p*E)V = W for V, Ã = tau*beta*A - 0.5*E 
%   or Ã = tau*beta*A - 0.5*E - UV^T if eqn.type == 'N'
%  Solves (Ã + p*E)^T*V = W for V, Ã = tau*beta*A - 0.5*E 
%   or Ã = tau*beta*A - 0.5*E - UV^T if eqn.type == 'T' 
%   (BDF scheme)
%
%
% Input:
%  eqn       structure containing equation data
%
%  opts      structure containing parameters for the algorithm
%
%  oper      contains function handles with operations for A and E
%
%  pc        contains shift parameter p
%
%  W         contains right hand side
%
% Output:
%  V         solution of the shifted system
%
%  eqn       structure containing equation data
%
%  opts      structure containing parameters for the algorithm
%
%  oper      contains function handles with operations for A and E

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

%% Check input

%% Initialize data
k = size(W, 2);
if eqn.haveUV
    m = size(eqn.U, 2);
end

%% Manipulate shift for BDF scheme
pc = (pc - 0.5) / (opts.bdf.tau * opts.bdf.beta);

%% preprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_pre(eqn, opts, oper);

%% solve shifted system
if eqn.haveUV %Perform Sherman-Morrison-Woodbury-trick
    V = oper.sol_ApE(eqn, opts,eqn.type,pc,eqn.type,...
        [W eqn.V] / (opts.bdf.tau * opts.bdf.beta),'N');
    SMW = V(:,k+1:end);
    V=V(:,1:k);
    V=V+SMW*((eye(m)-eqn.U'*SMW)\(eqn.U'*V));
else
    V = oper.sol_ApE(eqn, opts,eqn.type,pc,eqn.type,...
        W / (opts.bdf.tau * opts.bdf.beta),'N');
end

%% postprocess shifted solver
[eqn, opts, oper] = oper.sol_ApE_post(eqn, opts, oper);
