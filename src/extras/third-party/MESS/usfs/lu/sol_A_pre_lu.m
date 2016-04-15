function [ eqn, opts, oper ] = sol_A_pre_lu( eqn, opts, oper )
%% function pre initializes data and/or functions
%
% Input:
%    eqn    struct contains data for equations
%    
%    opts   struct contains parameters for the algorithm
%   
%    oper   struct contains function handles for operation with A
%
% Output:
% eqn
% opts
% oper

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
% Copyright (C) Jens Saak, Martin Koehler and others 
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
%

%% from lyapack
if ~isfield(eqn,'LP_L') || isempty(eqn.LP_L) %no previous LU-decomposition
[eqn.LP_L, eqn.LP_U, eqn.LP_a, eqn.LP_o, eqn.LP_S] = lu(eqn.A_,'vector');
end

