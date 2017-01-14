%
%% DEMO_NM_T
% Demo script for the 'T' (transposed) case of the algebraic Riccati
% equation with rank-k-update
%
%   (A - UV')'X + X(A - UV') - XBB'X + C'C = 0.
%
% This script shows the new data structures and use of an adapted Newton
% method (mess_lrnm).
%
%
% Steffen Werner, 2015-10-11.

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

%% Clear workspace.
clear;

%% Construction of data.
% The matrices are chosen as
%   A - a (n x n)-matrix, sparse
%   U - a (n x k)-matrix, dense
%   V - a (n x k)-matrix, dense
%   B - a (n x m)-matrix, dense
%   C - a (p x n)-matrix, dense

% Size of example data.
n = 100;
k = 2;
m = 3;
p = 1;

% Construction of all matrices.
A = sparse( 1:n, 1:n, -2 * ones( 1, n ), n, n ) ...
    + sparse( 2:n, 1:n-1, ones( 1, n - 1 ), n, n ) ...
    + sparse( 2:n, 1:n-1, ones( 1, n - 1 ), n, n )';

U = 1 / n * rand( n, k );
V = 1 / n * rand( n, k );

B = [ ones( n - m, m ); eye( m, m ) ];
C = ones( p, n );

%% Construction of MESS structure.
% Set data.
eqn.A_ = A;
eqn.B  = B;
eqn.C  = C;
eqn.U  = U;
eqn.V  = V;

% Set rank-k-update flag.
eqn.setUV  = 1;
eqn.UVsize = k;

% Set equation type. ('T' = transposed, 'N' = not transposed)
eqn.type = 'T';

% ADI settings.
opts.adi.maxiter         = 100;
opts.adi.restol          = 1.0e-10;
opts.adi.rctol           = 0;
opts.adi.info            = 1;
opts.adi.projection.freq = 0;
opts.adi.computeZ        = 1;

% ADI shifts.
opts.adi.shifts.l0 = 5;
opts.adi.shifts.kp = 10;
opts.adi.shifts.km = 5;
opts.adi.shifts.b0 = ones( n, 1 );

% NM settings.
opts.nm.maxiter = 50;
opts.nm.restol  = 1.0e-10;
opts.nm.rctol   = 1.0e-16;
opts.nm.info    = 1;
opts.nm.linesearch = 1;
opts.nm.accumulateRes = 1;
opts.nm.accumulateK = 1;
opts.nm.norm = 'fro';

% Set operations.
oper = operatormanager( 'default' );

%% Usage of the method.
[ Z, out ] = mess_lrnm( eqn, opts, oper );

%% Test of the solution.
abserr = norm( (A - U * V')' * (Z * Z') + (Z * Z') * (A - U * V') - (Z * Z') * (B * (B' * (Z * Z'))) + C' * C, 'fro' );
relerr = abserr / norm( Z * Z', 'fro' );
