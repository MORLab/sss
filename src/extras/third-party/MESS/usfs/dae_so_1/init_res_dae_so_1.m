function [ W, res0 ] = init_res_dae_so_1( eqn, opts, RHS)
%% function init_res initializes the low rank residual W and res0
%
%  Input:
%     eqn     structure contains data for A, B and C
%
%     opts    structure contains parameters for the algorithm
%
%     RHS     right hand side matrix
% 
%  Output:
%
%    W
%    res0
%
%   uses no other dae_so_1 function
%% check data in eqn structure

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
if (~isfield(eqn,'K_') || ~isnumeric(eqn.K_))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
end
if (~isfield(eqn,'isSym'))
    isSym = 0;
else
    isSym = eqn.isSym;
end
if ~isfield(eqn, 'nd')    || ~isnumeric(eqn.nd)
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end
if ~isfield(eqn,'type')
  eqn.type='N';
  warning('MESS:equation_type',['Unable to determine type of equation.'...
    'Falling back to type ''N''']);
end

rowK = size(eqn.K_, 1);
nd = eqn.nd;
na = rowK - nd;

%% check input Paramters
if (~isnumeric(RHS)) || (~ismatrix(RHS))
    error('MESS:error_arguments','RHS has to ba a matrix');
end
[rowG, colG] = size(RHS);


if(nd + na ~= rowG)
    error('MESS:error_arguments', 'number of rows of RHS'' has to coincide with dimension of M, D, K.');
end
%% compute low rank residual

if eqn.type == 'B'
    w = RHS(1 : nd, : ) - eqn.K_(1 : nd, nd + 1 : end) * ...
        (eqn.K_(nd + 1 : end, nd + 1 : end) \ (RHS(1 : nd, : )));
else
    if isSym
        w = RHS(1 : nd, : ) - eqn.K_(1 : nd, nd + 1 : end) * ...
            (eqn.K_(nd + 1 : end, nd + 1 : end) \ (RHS(1 : nd, : )));
    else
        w = RHS(1 : nd, : ) - eqn.K_(nd + 1 : end, 1 : nd)' * ...
            (eqn.K_(nd + 1 : end, nd + 1 : end)' \ (RHS(1 : nd, : )));
    end
end

W = [zeros(nd, colG); w];
%% compute res0
if opts.adi.LDL_T
    res0 = max(abs(eig(W' * W * eqn.S)));
else
    res0 = norm(W' * W, 2);
end

end

