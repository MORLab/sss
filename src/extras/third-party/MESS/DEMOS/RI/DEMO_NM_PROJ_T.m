%
%% DEMO_NM_PROJ_T
% Demo script for the 'T' (transposed) case of the algebraic Riccati
% equation with rank-k-update 
%
%   (A - UV')'X + X(A - UV') - XBB'X + C'C = 0.
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
% Loads Data.
[ E, A, B, C ] = getrail( 1 );

% Generates rank-k disturbances.
n = size( A, 1 );
k = 25;
U = 1 / n * rand( n, k );
V = 1 / n * rand( n, k );

%% Construction of MESS structure.
% Set data.
eqn.E_ = E;
eqn.A_ = A;
eqn.B  = B;
eqn.C  = C;
eqn.U  = U;
eqn.V  = V;

% Choice of Generalized or normal version.
eqn.haveE = 0;

% Set rank-k-update flag.
eqn.setUV  = 1;
eqn.UVsize = k;

% Set equation type. ('T' = transposed, 'N' = not transposed)
eqn.type = 'T';

% ADI settings.
opts.adi.maxiter         = 100;
opts.adi.restol          = 1.0e-14;
opts.adi.rctol           = 1.0e-16;
opts.adi.info            = 1;
opts.adi.projection.freq = 0;
opts.adi.computeZ        = 1;

% ADI shifts.
opts.adi.shifts.l0 = 25;
opts.adi.shifts.kp = 50;
opts.adi.shifts.km = 25;
opts.adi.shifts.b0 = ones( n, 1 );

% NM settings.
opts.nm.maxiter = 50;
opts.nm.restol  = 1.0e-08;
opts.nm.rctol   = 1.0e-16;
opts.nm.info    = 1;

% NM projection settings.
%opts.nm.projection      = 1;
opts.nm.projection.freq = 1;
opts.nm.res.maxiter     = 10;
opts.nm.res.tol         = 1.0e-06;
opts.nm.res.orth        = 1;
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
