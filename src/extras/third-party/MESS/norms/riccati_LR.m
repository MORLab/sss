function [ res ] = riccati_LR(W, DeltaK, opts, S, S_K )
%RICCATI_LR Norm of the Riccati residual
%   Using low rank formulation
%   R(X) = U * D * U^T
%   Since U^T * U * D is not symmetric ||U * D * U^T|| ~= ||U^T * U * D||
%   But the non-zero eigenvalues of (U * D * U^T) and (U^T * U * D) are
%   equal.

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

%% Compute norm

% 2-norm
if opts.nm.norm == 2
    if opts.adi.LDL_T
        if isempty(S_K)
            res = max(abs(eig([W DeltaK]'*[W * S -DeltaK])));
        else
            res = max(abs(eig([W DeltaK]'*[W * S -DeltaK * S_K])));
        end
    else
        res = max(abs(eig([W DeltaK]'*[W -DeltaK])));
    end
elseif strcmp(opts.nm.norm, 'fro')
    % Fromenius norm
    if opts.adi.LDL_T
        if isempty(S_K)
            res = norm(eig([W, DeltaK]' * [W * S, -DeltaK]), 'fro');
        else
            res = norm(eig([W, DeltaK]' * [W * S, -DeltaK * S_K]), 'fro');
        end
    else
        res = norm(eig([W, DeltaK]' * [W, -DeltaK]), 'fro');
    end
end

end

