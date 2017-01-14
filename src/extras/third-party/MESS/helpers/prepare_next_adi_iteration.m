function [ opts, out, stop ] = prepare_next_adi_iteration( opts, out, res, rc, outer_res, i)
% Evaluate stopping criteria of LRADI for exact and inexact case and check
% whether line search is necessary.

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

%% Evaluate stopping criteria
stop = 0;
if ~isempty(outer_res) && (outer_res(i) < opts.nm.restol)
    % Riccati tolerance reached, stop
    stop = 1;
elseif opts.adi.inexact
    if opts.adi.restol && (i > 2) && (res(i) > res(1))
        % Lyapunov residual is growing, stop and perform line search
        stop = 1;
        out.linesearch = 1;
    elseif opts.adi.restol && (res(i) < opts.adi.outer_tol) && ...
            ((res(i) >= 1e2 * opts.nm.restol) || ((i > 1) && ...
            (outer_res(i) >= outer_res(i - 1))))
        % Outer tolerance reached and not close to finish Newton iteration
        stop = 1;
    end
    % Outer tolerance not reached and Lyapunov residual is not growing or
    % outer tolerance is reached but Newton iteration is almost finished,
    % do NOT stop ADI iteration
elseif ((opts.adi.restol && (res(i) < opts.adi.restol)) || ...
        (opts.adi.rctol && (rc(i) < opts.adi.rctol))) && ...
        (isempty(outer_res) || ((res(i) >= 1e2 * opts.nm.restol) || ...
        ((i > 1) && (outer_res(i) >= outer_res(i - 1)))))
    % ADI tolarance reached, stop ADI iteration
    stop = 1;
    if ~isempty(outer_res) && outer_res(i) > outer_res(1)
        % Riccati residual is growing, perform line search
        out.linesearch = 1;
    end
elseif opts.adi.restol && (res(i) > res(1) * 1e2)
    % Lyapunov residual is growing and inexact ADI with line search probably
    % failed already, stop ADI iteration and restart Newton iteration with
    % exact ADI iteration
    stop = 1;
    out.restart = 1;
end
end

