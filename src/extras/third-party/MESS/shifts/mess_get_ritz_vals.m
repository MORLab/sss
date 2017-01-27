function [rw,  Hp, Hm, Vp, Vm, eqn, opts, oper] = mess_get_ritz_vals(eqn,opts,oper)
%% check data

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
if ~isfield(opts,'adi') || ~isstruct(opts.adi)
  error('MESS:control_data','ADI control structure opts.ADI missing.');
end
if ~isfield(opts.adi,'shifts') || ~isstruct(opts.adi.shifts)
    warning('MESS:control_data',['shift parameter control structure missing.', ...
        'Switching to default l0 = 25, kp = 50, km = 25.']);
    opts.adi.shifts.l0 = 25;
    opts.adi.shifts.kp = 50;
    opts.adi.shifts.km = 25;
else
    if ~isfield(opts.adi.shifts,'l0')||~isnumeric(opts.adi.shifts.l0)
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.adi.shifts.l0 field.', ...
            'Switching to default: 25']);
        opts.adi.shifts.l0 = 25;
    end
    if ~isfield(opts.adi.shifts,'kp')||~isnumeric(opts.adi.shifts.kp)
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.adi.shifts.kp field.', ...
            'Switching to default: 50']);
        opts.adi.shifts.kp = 50;
    end
    if ~isfield(opts.adi.shifts,'km')||~isnumeric(opts.adi.shifts.km)
        warning('MESS:control_data',...
            ['Missing or Corrupted opts.adi.shifts.km field.', ...
            'Switching to default: 25']);
        opts.adi.shifts.km = 25;
    end
end
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
[eqn, erg] = oper.init(eqn, opts, 'A','E');
if ~erg
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end

n = oper.size(eqn, opts);

if opts.adi.shifts.kp >= n, error('kp must be smaller than n!'); end
if opts.adi.shifts.km >= n, error('km must be smaller than n!'); end
if (2 * (opts.adi.shifts.l0) >= opts.adi.shifts.kp + opts.adi.shifts.km), ...
    error('2*l0 must be smaller than kp+km!'); end

if (~isfield(opts.adi.shifts, 'b0') || isempty(opts.adi.shifts.b0))
  opts.adi.shifts.b0 = ones(n, 1);
end


%% initialize data
opts.adi.shifts.b0 = (1 / norm(opts.adi.shifts.b0)) * opts.adi.shifts.b0;
rwp = [];
rwm = [];
rw = [];
Hp = [];
Vp = [];
Hm = [];
Vm = [];

%% estimate suboptimal ADI shift parameters
if opts.adi.shifts.kp > 0
  [Hp, Vp] = mess_arn(eqn, opts, oper, 'N');
  rwp = eig(Hp(1:opts.adi.shifts.kp, 1:opts.adi.shifts.kp));                 % =: R_+
  rw = [rw; rwp];
end

if opts.adi.shifts.km > 0
  [Hm, Vm] = mess_arn(eqn, opts, oper, 'I');
  rwm = ones(opts.adi.shifts.km, 1)./eig(Hm(1:opts.adi.shifts.km, ...
    1:opts.adi.shifts.km));     % =: 1 / R_-
  rw = [rw; rwm];                           % =: R
end
if any(real(rw) >= zeros(size(rw)))
  err_code = 1;
  %disp('These are the Ritz values computed by the Arnoldi process w.r.t. F:')
  %disp(rwp)
  %disp('These are the Ritz values computed by the Arnoldi process w.r.t. inv(F):')
  %disp(rwm)
%   disp(' ');
%   disp('####################################################################');
%   disp('WARNING in ''mess_para'': NON-STABLE RITZ VALUES DETECTED!!!')
%   disp(' ');
%   disp('This is quite a serious problem, that can be caused by  ');
%   disp('(i)   non-stable matrices F (Be sure that F is stable. ADI like');
%   disp('      methods only work for stable or antistable problems. If your');
%   disp('      Lyapunov equation is antistable, multiply it by -1.)');
%   disp('(ii)  matrices F that are stable but have an indefinite symmetric')
%   disp('      part (This is THE weak point of this algorithm. Try to work')
%   disp('      with the "reduced" Ritz values, i.e., the unstable values are')
%   disp('      simply removed. This is not an elegant measure but it may work.')
%   disp('      However, the convergence of ADI can be poor. This measure is')
%   disp('      taken automatically. Another measure might be to enlarge the')
%   disp('      values of kp or km, and run the program again.')
%   disp('(iii) matrices F with a negative definite, but ill-conditioned')
%   disp('      symmetric part (This is quite unlikely. The problem is')
%   disp('      caused by round-off errors).')
%   disp(' ')
%   disp('#####################################################################')
%   disp(' ');
%   disp(' ');
%   disp('NOTE: The unstable Ritz values will be ignored in the further computation!!! ');
%   disp(' ')
  %pause(3);

  rw  = rw(real(rw)<0);
end
end
