function [p, AV] = mess_projection_shifts(eqn, opts, oper, V, W, p_old)

%% Check data

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
if ~isfield(opts, 'adi') || ~isstruct(opts.adi)
    error('MESS:control_data', 'ADI control structure opts.ADI missing.');
end
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
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
[eqn, erg] = oper.init(eqn, opts, 'A', 'E');
if ~erg
    error('MESS:control_data', 'system data is not completely defined or corrupted');
end

L = length(p_old);
nV = size(V, 2);
nW = size(W, 2);
if L > 0 && any(p_old)
    if nV / nW ~= L
        error('MESS:control_data', 'V and W have inconsistend no. of columns');
    end
end

%% Initialize data
if L > 0 && any(p_old)
    T = zeros(L, L);
    K = zeros(1, L);
    D = [];
    Ir = eye(nW);
    iC = find(imag(p_old));
    iCh = iC(1 : 2 : end);
    iR = find(~imag(p_old));
    isubdiag = [iR; iCh];
    h = 1;
end

%% Process previous shifts
if L > 0 && any(p_old)
    while h <= L
        is = isubdiag(isubdiag < h);
        K(1, h) = 1;
        if isreal(p_old(h)) % real shift
            T(h, h) = p_old(h);
            if ~isempty(is)
                T(h, is) = 2 * p_old(h) * ones(1, length(is));
            end
            D = blkdiag(D, sqrt(-2 * p_old(h)));
            h = h + 1;
        else % complex conjugated pair of shifts
            rpc=real(p_old(h));
            ipc=imag(p_old(h));
            beta=rpc / ipc;
            T(h : h + 1, h : h + 1) = [3 * rpc, -ipc;
                                       ipc * (1 + 4 * beta^2), -rpc];
            if ~isempty(is)
                T(h : h+  1, is)=[4 * rpc; 
                                  4 * rpc * beta] * ones(1, length(is));
            end
            D = blkdiag(D, 2 * sqrt(-rpc) * [1,0; beta, sqrt(1 + beta^2)]);
            h = h + 2;
        end
    end
    S = kron(D \ (T * D), Ir); 
    K = kron(K * D, Ir); 
else
    S = 0;
    K = 1;
end

%% Compute projection matrices
[~, s, v] = svd(V' * V);
s = diag(s);
r = sum(s > eps * s(1) * nV);
st = v( : , 1 : r) * diag(1 ./ s(1 : r).^.5);
U = V * st;

%% Project V and compute Ritz values
if eqn.haveE
    E_V = oper.mul_E(eqn, opts, eqn.type, V, 'N');
    G = U' * E_V;
    H = U' * W * K * st + G * (S * st);
    G = G * st;
    p = eig(H, G);
%     AV = W * K + E_V * S;
else
    H = U' * (W * K) * st + U' * (V * S * st);
    p = eig(H);
%     AV = W * K + V * S;
end

%% Test
% if isfield(eqn, 'haveUV') && eqn.haveUV && (L > 0)
%     if eqn.type == 'T'
%         AV_ = oper.mul_A(eqn, opts, 'T', V, 'N') - eqn.V * (eqn.U' * V);
%     else
%         AV_ = oper.mul_A(eqn, opts, 'N', V, 'N') - eqn.U * (eqn.V' * V);
%     end
% else
%     AV_ = oper.mul_A(eqn, opts, eqn.type, V, 'N');
% end
% err = norm(AV - AV_);
% if ((err / norm(AV)) > 1e-10) && (err > 1e-10)
%     err
% end

%% Postprocess new shifts

% remove infinite values
p = p(isfinite(p));
% remove zeros
p = p(abs(p) > eps);
% make all shifts stable
p(real(p) > 0) = -p(real(p) > 0);
if ~isempty(p)
    % remove small imaginary pertubations
    small_imag = find(abs(imag(p)) ./ abs(p) < 1e-12);
    p(small_imag) = real(p(small_imag));
    % sort (s.t. compl. pairs are together)
    p=sort(p);
    l0 = min(opts.adi.shifts.l0, length(p));
%     l0 = min(min(nW, opts.adi.shifts.l0), length(p));
    p = mess_mnmx(p, l0);
end
