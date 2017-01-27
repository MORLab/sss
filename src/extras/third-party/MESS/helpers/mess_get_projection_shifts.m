function [ opts, l ] = mess_get_projection_shifts( eqn, opts, oper, Z, W, D)
% compute next projection shifts and update shift vectors if shift
% computation method is 'projection', if method is not 'projection' nothig
% is done

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
k = size(eqn.G, 2);
if ~isfield(opts.adi,'shifts') || ~isstruct(opts.adi.shifts)
    warning('MESS:control_data',['shift parameter control structure missing.', ...
        'Switching to default l0 = 25.']);
    opts.adi.shifts.l0 = 25;
else
    if ~isfield(opts.adi.shifts,'l0')||~isnumeric(opts.adi.shifts.l0)
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.adi.shifts.l0 field.', ...
            'Switching to default: 25']);
        opts.adi.shifts.l0 = 25;
    end
end

%% Compute new shifts
if isfield(opts.adi.shifts, 'method') && ...
        strcmp(opts.adi.shifts.method, 'projection')
    if ~isfield(opts.adi.shifts, 'used_shifts') || ...
            isempty(opts.adi.shifts.used_shifts)
        opts.adi.shifts.used_shifts = opts.adi.shifts.p;
    else
        opts.adi.shifts.used_shifts = ...
            [opts.adi.shifts.used_shifts; opts.adi.shifts.p];
        if (size(opts.adi.shifts.used_shifts, 1) > opts.adi.shifts.l0) && ...
                imag(opts.adi.shifts.used_shifts(end - opts.adi.shifts.l0 + 1)) && ...
                (abs(opts.adi.shifts.used_shifts(end - opts.adi.shifts.l0 + 1) ...
                -conj(opts.adi.shifts.used_shifts(end - opts.adi.shifts.l0))) < eps)
            % don't cut between pair of complex shifts
            opts.adi.shifts.used_shifts = ...
                opts.adi.shifts.used_shifts(end - opts.adi.shifts.l0 : end);
        elseif (size(opts.adi.shifts.used_shifts, 1) > opts.adi.shifts.l0)
            opts.adi.shifts.used_shifts = ...
                opts.adi.shifts.used_shifts(end - opts.adi.shifts.l0 + 1 : end);
        end
    end
    if isfield(oper,'get_ritz_vals')
        if opts.adi.LDL_T
            % scale columns of Z (L) as in original non LDL^T formulation
            len = size(opts.adi.shifts.used_shifts, 1) * k - 1;
            p = oper.get_ritz_vals(eqn, opts, oper, ...
                Z( : ,end - len : end) * ...
                kron(diag(sqrt(D(end - size(opts.adi.shifts.used_shifts, 1) + 1: end))), ...
                eye(k)), W, opts.adi.shifts.used_shifts);
        else
            p = oper.get_ritz_vals(eqn, opts, oper, ...
                Z( : ,end - (size(opts.adi.shifts.used_shifts, 1) * k) + 1 : end), ...
                W, opts.adi.shifts.used_shifts);
        end
    else
        if opts.adi.LDL_T
            % scale columns of Z (L) as in original non LDL^T formulation
            len = size(opts.adi.shifts.used_shifts, 1) * k - 1;
            p = mess_projection_shifts(eqn, opts, oper, ...
                Z( : ,end - len : end) * ...
                kron(diag(sqrt(D(end - size(opts.adi.shifts.used_shifts, 1) + 1: end))), ...
                eye(k)), W, opts.adi.shifts.used_shifts);
        else
            p = mess_projection_shifts(eqn, opts, oper, ...
                Z( : ,end - (size(opts.adi.shifts.used_shifts, 1) * k) + 1 : end), ...
                W, opts.adi.shifts.used_shifts);
        end
    end
    if ~isempty(p)
        opts.adi.shifts.p = p;
        % else could not compute new shifts, use previous ones
        % again
    end
end
l = length(opts.adi.shifts.p);
end

