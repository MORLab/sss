function n = size_dae3_so(eqn, opts, oper)
% function n = size_dae3_so(eqn, opts)
%
% This function returns the number of rows of the implicitly projected A 
% matrix of the index-3 system.
%
%    Inputs:
%
%    eqn     structure containing the system data
%
%    Output:
%
%    n       size of the implicitly projected A matrix

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

n= 2*size(eqn.M_,1);
end
