function n = size_dae_so_1(eqn, ~)

% function n = size_dae_so_1(eqn, opts)
%
% This function returns the number of rows of matrix A in the Schur
% complement formulation of the index-1 system.
%
%    Inputs:
%
%    eqn     structure containing the system data
%
%    Output:
%
%    n       size of the implicitly projected A matrix
%
% This function does not use other dae_1_so_1 functions.

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
if ~isfield(eqn, 'nd')    || ~isnumeric(eqn.nd)
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end
n = 2 * eqn.nd;
end
