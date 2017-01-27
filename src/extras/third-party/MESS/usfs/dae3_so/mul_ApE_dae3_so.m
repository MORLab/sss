function C = mul_ApE_dae3_so(eqn, opts, opA,p,opE, B, opB)
%% function mul_ApE_dae3_so computes C = (opA(A_) + p*opE(E_))*opB(B)
%
%
% Input:
%   eqn         structure contains matrices
%
%   eqn.typeE   specifies whether E_ is Identity or not
%               typeE = 0 E_ is Identity
%               typeE = 1 E_ is not Identity
%
%   opts        struct contains parameters for the algorithm
%
%   opA         character specifies the form of opA(A_)
%               opA = 'N' performs (A_+pc*opE(E_))*opB(B)
%               opA = 'T' performs (A_'+pc*opE(E_))*opB(B)
%
%   p           scalar value
%
%   opE         character specifies the form of opE(E_)
%               opE = 'N' performs (opA(A_)+pc*E_)*opB(B)
%               opE = 'T' performs (opA(A_)+pc*E_')*opB(B)
%
%   B           m-x-p matrix
%
%   opB         character specifies the form of opB(B)
%               opB = 'N' performs (opA(A_)+pc*opE(E_))*B
%               opB = 'T' performs (opA(A_)+pc*opE(E_))*B'
%
% Output:
% B = (opA(A_)+pc*opE(E_))*opB(B)
%
%   uses no other dae3_so function
%

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

%% check input Paramters
if (~ischar(opA) || ~ischar(opE) || ~ischar(opB))
    error('MESS:error_arguments', 'opA, opE or opB is not a char');
end

opA = upper(opA); opE = upper(opE); opB = upper(opB);

if(~(opA == 'N' || opA == 'T'))
    error('MESS:error_arguments', 'opA is not ''N'' or ''T''');
end

if(~(opE == 'N' || opE == 'T'))
    error('MESS:error_arguments', 'opE is not ''N'' or ''T''');
end

if(~(opB == 'N' || opB == 'T'))
    error('MESS:error_arguments', 'opB is not ''N'' or ''T''');
end

if(~isnumeric(p))
   error('MESS:error_arguments','p is not numeric'); 
end

if (~isnumeric(B)) || (~ismatrix(B))
    error('MESS:error_arguments','B has to ba a matrix');
end

%% check data in eqn structure
for mat='MDKG'
    if(~isfield(eqn, sprintf('%c_',mat)) || ~eval(sprintf('isnumeric(eqn.%c_)',mat)))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

nv = size(eqn.M_,1);
np = size(eqn.G_,1);

[rowB,colB] = size(B);

if(opB == 'N')
    if (rowB ~= 2*nv + np)
        B = [B; zeros(2 * nv + np - rowB, colB)];
    end
else
    if (colB ~= 2*nv + np)
        B = [B, zeros(rowB, 2 * nv + np - colB)];
    end
end

%% compute C = (A + p * E) * B
if eqn.haveE
    %% perfom solve operations for E ~= Identity
    switch opA
        case 'N'
            switch opE
                case 'N'
                    switch opB
                        case 'N'
                            C1 = p * B(1 : nv, :) + B(nv + 1 : 2 * nv, :);
                            C2 = eqn.K_ * B(1 : nv, :) ...
                                + (eqn.D_ + p * eqn.M_) * B(nv + 1 : 2 * nv, :) ...
                                + eqn.G_' * B(2 * nv + 1 : end, :);
                            C3 = eqn.G_ * B(1 : nv, :);
                            C = [C1; C2; C3];
                        case 'T'
                            C1 = p * B( : , 1 : nv)' + B( : , nv + 1 : 2 * nv)';
                            C2 = eqn.K_ * B( : , 1 : nv)' ...
                                + (eqn.D_ + p * eqn.M_) * B( : , nv + 1 : 2 * nv)' ...
                                + eqn.G_' * B( : , 2 * nv + 1 : end)';
                            C3 = eqn.G_ * B( : , 1 : nv)';
                            C = [C1; C2; C3];
                    end
                case 'T'
                    switch opB
                        case 'N'
                            C1 = p * B(1 : nv, :) + B(nv + 1 : 2 * nv, :);
                            C2 = eqn.K_ * B(1 : nv, :) ...
                                + (eqn.D_ + p * eqn.M_') * B(nv + 1 : 2 * nv, :) ...
                                + eqn.G_' * B(2 * nv + 1 : end, :);
                            C3 = eqn.G_ * B(1 : nv, :);
                            C = [C1; C2; C3];
                        case 'T'
                            C1 = p * B( : , 1 : nv)' + B( : , nv + 1 : 2 * nv)';
                            C2 = eqn.K_ * B( : , 1 : nv)' ...
                                + (eqn.D_ + p * eqn.M_') * B( : , nv + 1 : 2 * nv)' ...
                                + eqn.G_' * B( : , 2 * nv + 1 : end)';
                            C3 = eqn.G_ * B( : , 1 : nv)';
                            C = [C1; C2; C3];
                    end
            end
        case 'T'
            switch opE
                case 'N'
                    switch opB
                        case 'N'
                            C1 = p * B(1 : nv, :) + eqn.K_' * B(nv + 1 : 2 * nv, :) ...
                                + eqn.G_' * B(2 * nv + 1 : end, :);
                            C2 = B(1 : nv, :) ...
                                + (eqn.D_' + p * eqn.M_) * B(nv + 1 : 2 * nv, :);
                            C3 = eqn.G_ * B(nv + 1 : 2 * nv, :);
                            C = [C1; C2; C3];
                        case 'T'
                            C1 = p * B( : , 1 : nv)' + eqn.K_' * B( : , nv + 1 : 2 * nv)' ...
                                + eqn.G_' * B( : , 2 * nv + 1 : end)';
                            C2 = B( : , 1 : nv)' ...
                                + (eqn.D_' + p * eqn.M_) * B( : , nv + 1 : 2 * nv)';
                            C3 = eqn.G_ * B( : , nv + 1 : 2 * nv)';
                            C = [C1; C2; C3];
                    end
                case 'T'
                    switch opB
                        case 'N'
                            C1 = p * B(1 : nv, :) + eqn.K_' * B(nv + 1 : 2 * nv, :) ...
                                + eqn.G_' * B(2 * nv + 1 : end, :);
                            C2 = B(1 : nv, :) ...
                                + (eqn.D_' + p * eqn.M_') * B(nv + 1 : 2 * nv, :);
                            C3 = eqn.G_ * B(nv + 1 : 2 * nv, :);
                            C = [C1; C2; C3];
                        case 'T'
                            C1 = p * B( : , 1 : nv)' + eqn.K_' * B( : , nv + 1 : 2 * nv)' ...
                                + eqn.G_' * B( : , 2 * nv + 1 : end)';
                            C2 = B( : , 1 : nv)' ...
                                + (eqn.D_' + p * eqn.M_') * B( : , nv + 1 : 2 * nv)';
                            C3 = eqn.G_ * B( : , nv + 1 : 2 * nv)';
                            C = [C1; C2; C3];
                    end
            end
    end
else
    %% perfom solve operations for E = Identity
    switch opA
        case 'N'
            switch opB
                case 'N'
                    C1 = p * B(1 : nv, :) + B(nv + 1 : 2 * nv, :);
                    C2 = eqn.K_ * B(1 : nv, :) ...
                        + (eqn.D_ + p * speye(nv, nv)) * B(nv + 1 : 2 * nv, :) ...
                        + eqn.G_' * B(2 * nv + 1 : end, :);
                    C3 = eqn.G_ * B(1 : nv, :);
                    C = [C1; C2; C3];
                case 'T'
                    C1 = p * B( : , 1 : nv)' + B( : , nv + 1 : 2 * nv)';
                    C2 = eqn.K_ * B( : , 1 : nv)' ...
                        + (eqn.D_ + p * speye(nv, nv)) * B( : , nv + 1 : 2 * nv)' ...
                        + eqn.G_' * B( : , 2 * nv + 1 : end)';
                    C3 = eqn.G_ * B( : , 1 : nv)';
                    C = [C1; C2; C3];
            end
        case 'T'
            switch opB
                case 'N'
                    C1 = p * B(1 : nv, :) + eqn.K_' * B(nv + 1 : 2 * nv, :) ...
                        + eqn.G_' * B(2 * nv + 1 : end, :);
                    C2 = B(1 : nv, :) ...
                        + (eqn.D_' + p * speye(nv, nv)) * B(nv + 1 : 2 * nv, :);
                    C3 = eqn.G_ * B(nv + 1 : 2 * nv, :);
                    C = [C1; C2; C3];
                case 'T'
                    C1 = p * B( : , 1 : nv)' + eqn.K_' * B( : , nv + 1 : 2 * nv)' ...
                        + eqn.G_' * B( : , 2 * nv + 1 : end)';
                    C2 = B( : , 1 : nv)' ...
                        + (eqn.D_' + p * speye(nv, nv)) * B( : , nv + 1 : 2 * nv)';
                    C3 = eqn.G_ * B( : , nv + 1 : 2 * nv)';
                    C = [C1; C2; C3];
            end
    end
end
if opB == 'N'
    C = C(1 : rowB, : );
else
    C = C(1 : colB, : );
end
end
