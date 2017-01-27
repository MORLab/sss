function Z=mess_column_compression(Z,opts)
%          Computes a compressed representation of Z using the rank
%          revealing SVD.
%
%   Input
%       Z             Matrix of interest
%                 
%       opts          structure with folowing fields
%       opts.ccTOL    the truncation tolerance used in the RRQRprocess.
%                     Note that in the context of factorized solution of matrix
%                     equations it is sufficient to choose
%                     sqrt(eps) here since in the product this will
%                     lead to a truncation error eps. Only do this
%                     when the product is used in subsequent
%                     computations, though.
%
%   Output
%       Z             compressed low rank factor
%
% author  Jens Saak
% date    2012/07/01

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

if ~isfield(opts,'comprType')||(isempty(opts.comprType))
  opts.comprType=0; 
end
if ~isfield(opts,'ccTOL')||(isempty(opts.ccTOL))
  opts.ccTOL=eps; 
end

if(issparse(Z)) 
  % This is just a safety measure that is hopefully never executed
  Z=full(Z);
  warning('MESS:dense',['Converting low rank factor to dense format. ' ...
                      'This should never be necessary.']')
end

[~,D,P] = svd(Z');

L = D';
D = diag(D);
l = length(find(D > opt.ccTOL*D(1)));
Z = P * L(:,1:l);
