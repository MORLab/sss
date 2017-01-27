function [G, S]=mess_LDL_column_compression(G, opG, S, opS, tol)
%   Computes a compressed representation of G * S * G' using the rank
%          revealing SVD.
% author  Björn Baran
% date    2015/08/24
%
%   Input
%       G,S         Matrices of interest
%       opG, opS    Character specifing whether G or S should be
%                   transposed
%       G,S         Matrices of interest
%
%   Output
%       tol         Tolerance
%       G, S        compressed low rank factors

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


%% Compute QR decomposition
if opG == 'N'
    [Q, R, e] = qr(G, 0);
else
    [Q, R, e] = qr(G', 0);
end
e(e) = 1 : length(e);
R = R(:,e);

%% Compute eigendecomposition

if opS == 'N'
    RSR = R * S * R';
else
    RSR = R * S' * R';
end
RSR = (RSR + RSR') / 2;
[V, D] = eig(RSR);
D = diag(D);

%% Compute rank and truncate
r = abs(D) > tol * max(abs(D));
V = V(  : , r);
if opG == 'N'
    G = Q * V;
else
    G = (Q * V)';
end
S = diag(D(r));
