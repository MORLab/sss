%
%% DEMO_NM_GE_T
% Demo script for the 'T' (transposed) case of the generalized algebraic 
% Riccati equation with rank-k-update 
%
%   (A - UV')'XE + E'X(A - UV') - E'XBB'XE + C'C = 0.
%
% using projection methods.
% This script shows the new data structures and use of an adapted Newton
% method (mess_lrnm).
%
%
% Steffen Werner, 2015-10-10.

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
% Size of example data.
n = 100;
k = 2;
m = 3;
p = 1;

% Construction of all matrices.
A = sparse( 1:n, 1:n, -2 * ones( 1, n ), n, n ) ...
    + sparse( 2:n, 1:n-1, ones( 1, n - 1 ), n, n ) ...
    + sparse( 2:n, 1:n-1, ones( 1, n - 1 ), n, n )';

E = sparse( 1:n, 1:n, 3 * ones( 1, n ), n, n )...
    + sparse( 2:n, 1:n-1, ones( 1, n - 1 ), n, n ) ...
    + sparse( 2:n, 1:n-1, ones( 1, n - 1 ), n, n )';

U = 1 / n * rand( n, k );
V = 1 / n * rand( n, k );

B = [ ones( n - m, m ); eye( m, m ) ];
C = ones( p, n );

%% Construction of MESS structure.
% Set data.
eqn.E_ = E;
eqn.A_ = A;
eqn.B  = B;
eqn.C  = C;
eqn.U  = U;
eqn.V  = V;

% Choice of Generalized or normal version.
eqn.haveE = 1;

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
opts.nm.restol  = 1.0e-08;
opts.nm.rctol   = 1.0e-16;
opts.nm.info    = 1;
opts.nm.norm = 'fro';

opts.nm.linesearch = 1;
opts.nm.accumulateRes = 1;
opts.nm.accumulateK = 1;

% Set operations.
oper = operatormanager( 'default' );

%% Usage of the method.
[ Z, out ] = mess_lrnm( eqn, opts, oper );

%% Test of the solution.
abserr = norm( (A - U * V')' * (Z * Z') * E + E' * (Z * Z') * (A - U * V') - E' * (Z * Z') * (B * (B' * (Z * Z'))) * E + C' * C, 'fro' );
relerr = abserr / norm( (E' * Z) * Z', 'fro' );
