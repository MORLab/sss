%
%% DEMO_RI_GE_T2
% Demo script for the 'T' (transposed) case of the H-infinity algebraic 
% Riccati equation
%
%   A'XE + E'XA + E'X(B1B1' - B2B2')XE + C'C = 0.
%
% Demonstrates the low-rank Riccati iteration with M.E.S.S. funtions and
% structures.
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
% Loads Data.
[ E, A, B, C ] = getrail( 1 );

B1 = B(:,1:2);
B2 = B(:,3:end);

n = size( A, 1 );

%% Construction of MESS structure.
% Set data.
eqn.E_ = E;
eqn.A_ = A;
eqn.B1 = B1;
eqn.B2 = B2;
eqn.C  = C;

% Choice of Generalized or normal version.
eqn.haveE = 1;

% Set equation type. ('T' = transposed, 'N' = not transposed)
eqn.type = 'T';

% ADI settings.
opts.adi.maxiter         = 200;
opts.adi.restol          = 1.0e-10;
opts.adi.rctol           = 0;
opts.adi.info            = 0;
opts.adi.projection.freq = 0;
opts.adi.computeZ        = 1;
opts.adi.accumulateK     = 0;
opts.adi.accumulateDeltaK     = 0;

% ADI shifts.
opts.adi.shifts.l0 = 5;
opts.adi.shifts.kp = 10;
opts.adi.shifts.km = 5;
opts.adi.shifts.b0 = ones( n, 1 );

% NM settings.
opts.nm.maxiter = 50;
opts.nm.restol  = 1.0e-10;
opts.nm.rctol   = 1.0e-12;
opts.nm.info    = 1;

% NM projection settings.
opts.nm.projection      = [];
opts.nm.projection.freq = 0;
opts.nm.res.maxiter     = 10;
opts.nm.res.tol         = 1.0e-06;
opts.nm.res.orth        = 1;

opts.nm.linesearch = 0;
opts.nm.norm = 'fro';
opts.nm.accumulateRes = 1;

opts.nm.inexact = 'linear';
opts.nm.tau = 0.1;

% RI settings.
opts.ri.maxiter = 2;
opts.ri.restol  = 1.0e-08;
opts.ri.colctol = 1.0e-10;
opts.ri.info    = 1;

% Set operations.
oper = operatormanager( 'default' );

%% Usage of the method.
[ Z, out ] = mess_lrri( eqn, opts, oper );

%% Test of the solution.
abserr = norm( A' * (Z * Z') * E + E' * (Z * Z') * A + E' * (Z * Z') * (B1 * B1' - B2 * B2') * (Z * Z') * E + C' * C, 2 );
relerr = abserr / norm( (E' * Z) * Z', 2 );
