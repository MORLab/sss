function C = mul_ApE_dae_so_1(eqn, opts, opA, p, opE, B, opB)

%% function mul_A mul_ApE_dae_so_1 operation C = (opA(A_)+pc*opE(E_))*opB(B)
%
% Input:
%   eqn     structure contains matrices
%
%   opts    struct contains parameters for the algorithm
%
%   opA     character specifies the form of opA(A_)
%           opA = 'N' performs (A_+pc*opE(E_))*opB(B)
%           opA = 'T' performs (A_'+pc*opE(E_))*opB(B)
%
%   p       scalar value
%
%   opE     character specifies the form of opE(E_)
%           opE = 'N' performs (opA(A_)+pc*E_)*opB(B)
%           opE = 'T' performs (opA(A_)+pc*E_')*opB(B)
%
%   B       m-x-p matrix
%
%   opB     character specifies the form of opB(B)
%           opB = 'N' performs (opA(A_)+pc*opE(E_))*B
%           opB = 'T' performs (opA(A_)+pc*opE(E_))*B'
%
% Output:
% B = (opA(A_)+pc*opE(E_))*opB(B)
%
%   uses no other dae_so_1 function

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
if (~isfield(eqn,'K_') || ~isnumeric(eqn.K_))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
end
if (~isfield(eqn,'M_') || ~isnumeric(eqn.M_))
    error('MESS:equation_data',...
        'Empty or Corrupted field M detected in equation structure.')
elseif (~isfield(eqn,'D_') || ~isnumeric(eqn.D_))
    error('MESS:equation_data',...
        'Empty or Corrupted field D detected in equation structure.')
end
if (~isfield(eqn,'isSym'))
    isSym = 0;
else
    isSym = eqn.isSym;
end
if ~isfield(eqn, 'nd')    || ~isnumeric(eqn.nd)
    error('MESS:nd',...
        'Missing or Corrupted nd field detected in equation structure.');
end
if ~isfield(eqn,'haveE'), eqn.haveE=0; end

nd = eqn.nd;

if(opB == 'N')
    rowB = size(B, 1);
else
    rowB = size(B, 2);
end

if(2 * nd ~= rowB)
    error('MESS:error_arguments','Rows of A differs from rows of B');
end

%% compute C = (A + p * E) * B
if isSym
    %% Assume symmetric data
    if eqn.haveE
        %% perfom solve operations for E ~= Identity
        switch opE
            case 'N'
                switch opB
                    case 'N'
                        C1 = B(1 : nd, :) + p * B(nd + 1 : end, :);
                        C2 = p * (eqn.M_(1 : nd, 1 : nd) * B(1 : nd, :)...
                            + eqn.D_(1 : nd, 1 : nd) * B(nd + 1 : end, :)) ...
                            - eqn.K_(1 : nd, 1 : nd) * B(nd + 1 : end, :) ...
                            + eqn.K_(1 : nd, nd + 1 : end) * ...
                            (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                            \ (eqn.K_(nd + 1 : end, 1 : nd) * B(nd + 1 : end, :)));
                        C = [C1; C2];
                    case 'T'
                        C1 = B( : ,1 : nd)' + p * B( : ,nd + 1 : end)';
                        C2 = p * (eqn.M_(1 : nd, 1 : nd) * B( : ,1 : nd)'...
                            + eqn.D_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)') ...
                            - eqn.K_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)' ...
                            + eqn.K_(1 : nd, nd + 1 : end) * ...
                            (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                            \ (eqn.K_(nd + 1 : end, 1 : nd) * B( : ,nd + 1 : end)'));
                        C = [C1; C2];
                end
            case 'T'
                switch opB
                    case 'N'
                        C1 = B(1 : nd, :) + p * eqn.M_(1 : nd, 1 : nd) * B(nd + 1 : end, :);
                        C2 = p * (B(1 : nd, :)...
                            + eqn.D_(1 : nd, 1 : nd) * B(nd + 1 : end, :)) ...
                            - eqn.K_(1 : nd, 1 : nd) * B(nd + 1 : end, :) ...
                            + eqn.K_(1 : nd, nd + 1 : end) * ...
                            (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                            \ (eqn.K_(nd + 1 : end, 1 : nd) * B(nd + 1 : end, :)));
                        C = [C1; C2];
                    case 'T'
                        C1 = B( : ,1 : nd)' + p * eqn.M_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)';
                        C2 = p * (B( : ,1 : nd)'...
                            + eqn.D_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)') ...
                            - eqn.K_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)' ...
                            + eqn.K_(1 : nd, nd + 1 : end) * ...
                            (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                            \ (eqn.K_(nd + 1 : end, 1 : nd) * B( : ,nd + 1 : end)'));
                        C = [C1; C2];
                end
        end
    else
        %% perfom solve operations for E = Identity
        switch opB
            case 'N'
                C1 = (p + 1) * B(1 : nd, :);
                C2 = p * B(nd + 1 : end, :) ...
                    - eqn.K_(1 : nd, 1 : nd) * B(nd + 1 : end, :) ...
                    + eqn.K_(1 : nd, nd + 1 : end) * ...
                    (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                    \ (eqn.K_(nd + 1 : end, 1 : nd) * B(nd + 1 : end, :)));
                C = [C1; C2];
            case 'T'
                C1 = ( p + 1) * B( : ,1 : nd)';
                C2 = p * B( : ,nd + 1 : end)' ...
                    - eqn.K_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)' ...
                    + eqn.K_(1 : nd, nd + 1 : end) * ...
                    (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                    \ (eqn.K_(nd + 1 : end, 1 : nd) * B( : ,nd + 1 : end)'));
                C = [C1; C2];
        end
    end
else
    %% Do not assume symmetric data
    if eqn.haveE
        %% perfom solve operations for E ~= Identity
        switch opA
            case 'N'
                switch opE
                    case 'N'
                        switch opB
                            case 'N'
                                C1 = B(1 : nd, :) + p * B(nd + 1 : end, :);
                                C2 = p * (eqn.M_(1 : nd, 1 : nd) * B(1 : nd, :)...
                                    + eqn.D_(1 : nd, 1 : nd) * B(nd + 1 : end, :)) ...
                                    - eqn.K_(1 : nd, 1 : nd) * B(nd + 1 : end, :) ...
                                    + eqn.K_(1 : nd, nd + 1 : end) * ...
                                    (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                                    \ (eqn.K_(nd + 1 : end, 1 : nd) * B(nd + 1 : end, :)));
                                C = [C1; C2];
                            case 'T'
                                C1 = B( : ,1 : nd)' + p * B( : ,nd + 1 : end)';
                                C2 = p * (eqn.M_(1 : nd, 1 : nd) * B( : ,1 : nd)'...
                                    + eqn.D_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)') ...
                                    - eqn.K_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)' ...
                                    + eqn.K_(1 : nd, nd + 1 : end) * ...
                                    (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                                    \ (eqn.K_(nd + 1 : end, 1 : nd) * B( : ,nd + 1 : end)'));
                                C = [C1; C2];
                        end
                    case 'T'
                        switch opB
                            case 'N'
                                C1 = B(1 : nd, :) + p * eqn.M_(1 : nd, 1 : nd)' * B(nd + 1 : end, :);
                                C2 = p * (B(1 : nd, :)...
                                    + eqn.D_(1 : nd, 1 : nd)' * B(nd + 1 : end, :)) ...
                                    - eqn.K_(1 : nd, 1 : nd) * B(nd + 1 : end, :) ...
                                    + eqn.K_(1 : nd, nd + 1 : end) * ...
                                    (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                                    \ (eqn.K_(nd + 1 : end, 1 : nd) * B(nd + 1 : end, :)));
                                C = [C1; C2];
                            case 'T'
                                C1 = B( : ,1 : nd)' + p * eqn.M_(1 : nd, 1 : nd)' * B( : ,nd + 1 : end)';
                                C2 = p * (B( : ,1 : nd)'...
                                    + eqn.D_(1 : nd, 1 : nd)' * B( : ,nd + 1 : end)') ...
                                    - eqn.K_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)' ...
                                    + eqn.K_(1 : nd, nd + 1 : end) * ...
                                    (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                                    \ (eqn.K_(nd + 1 : end, 1 : nd) * B( : ,nd + 1 : end)'));
                                C = [C1; C2];
                        end
                end
            case 'T'
                switch opE
                    case 'N'
                        switch opB
                            case 'N'
                                C1 = B(1 : nd, :) + p * B(nd + 1 : end, :);
                                C2 = p * (eqn.M_(1 : nd, 1 : nd) * B(1 : nd, :)...
                                    + eqn.D_(1 : nd, 1 : nd) * B(nd + 1 : end, :)) ...
                                    - eqn.K_(1 : nd, 1 : nd)' * B(nd + 1 : end, :) ...
                                    + eqn.K_(nd + 1 : end, 1 : nd)' * ...
                                    (eqn.K_(nd + 1 : end, nd + 1 : end)' ...
                                    \ (eqn.K_(1 : nd, nd + 1 : end)' * B(nd + 1 : end, :)));
                                C = [C1; C2];
                            case 'T'
                                C1 = B( : ,1 : nd)' + p * B( : ,nd + 1 : end)';
                                C2 = p * (eqn.M_(1 : nd, 1 : nd) * B( : ,1 : nd)'...
                                    + eqn.D_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)') ...
                                    - eqn.K_(1 : nd, 1 : nd)' * B( : ,nd + 1 : end)' ...
                                    + eqn.K_(nd + 1 : end, 1 : nd)' * ...
                                    (eqn.K_(nd + 1 : end, nd + 1 : end)' ...
                                    \ (eqn.K_(1 : nd, nd + 1 : end)' * B( : ,nd + 1 : end)'));
                                C = [C1; C2];
                        end
                    case 'T'
                        switch opB
                            case 'N'
                                C1 = B(1 : nd, :) + p * eqn.M_(1 : nd, 1 : nd)' * B(nd + 1 : end, :);
                                C2 = p * (B(1 : nd, :)...
                                    + eqn.D_(1 : nd, 1 : nd)' * B(nd + 1 : end, :)) ...
                                    - eqn.K_(1 : nd, 1 : nd)' * B(nd + 1 : end, :) ...
                                    + eqn.K_(nd + 1 : end, 1 : nd)' * ...
                                    (eqn.K_(nd + 1 : end, nd + 1 : end)' ...
                                    \ (eqn.K_(1 : nd, nd + 1 : end)' * B(nd + 1 : end, :)));
                                C = [C1; C2];
                            case 'T'
                                C1 = B( : ,1 : nd)' + p * eqn.M_(1 : nd, 1 : nd)' * B( : ,nd + 1 : end)';
                                C2 = p * (B( : ,1 : nd)'...
                                    + eqn.D_(1 : nd, 1 : nd)' * B( : ,nd + 1 : end)') ...
                                    - eqn.K_(1 : nd, 1 : nd)' * B( : ,nd + 1 : end)' ...
                                    + eqn.K_(nd + 1 : end, 1 : nd)' * ...
                                    (eqn.K_(nd + 1 : end, nd + 1 : end)' ...
                                    \ (eqn.K_(1 : nd, nd + 1 : end)' * B( : ,nd + 1 : end)'));
                                C = [C1; C2];
                        end
                end
        end
    else
        %% perfom solve operations for E = Identity
        switch opA
            case 'N'
                switch opB
                    case 'N'
                        C1 = (p + 1) * B(1 : nd, :);
                        C2 = p * B(nd + 1 : end, :) ...
                            - eqn.K_(1 : nd, 1 : nd) * B(nd + 1 : end, :) ...
                            + eqn.K_(1 : nd, nd + 1 : end) * ...
                            (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                            \ (eqn.K_(nd + 1 : end, 1 : nd) * B(nd + 1 : end, :)));
                        C = [C1; C2];
                    case 'T'
                        C1 = ( p + 1) * B( : ,1 : nd)';
                        C2 = p * B( : ,nd + 1 : end)' ...
                            - eqn.K_(1 : nd, 1 : nd) * B( : ,nd + 1 : end)' ...
                            + eqn.K_(1 : nd, nd + 1 : end) * ...
                            (eqn.K_(nd + 1 : end, nd + 1 : end) ...
                            \ (eqn.K_(nd + 1 : end, 1 : nd) * B( : ,nd + 1 : end)'));
                        C = [C1; C2];
                end
            case 'T'
                switch opB
                    case 'N'
                        C1 = (p + 1) * B(1 : nd, :);
                        C2 = p * B(nd + 1 : end, :) ...
                            - eqn.K_(1 : nd, 1 : nd)' * B(nd + 1 : end, :) ...
                            + eqn.K_(nd + 1 : end, 1 : nd)' * ...
                            (eqn.K_(nd + 1 : end, nd + 1 : end)' ...
                            \ (eqn.K_(1 : nd, nd + 1 : end)' * B(nd + 1 : end, :)));
                        C = [C1; C2];
                    case 'T'
                        C1 = (p + 1) * B( : ,1 : nd)';
                        C2 = p * B( : ,nd + 1 : end)' ...
                            - eqn.K_(1 : nd, 1 : nd)' * B( : ,nd + 1 : end)' ...
                            + eqn.K_(nd + 1 : end, 1 : nd)' * ...
                            (eqn.K_(nd + 1 : end, nd + 1 : end)' ...
                            \ (eqn.K_(1 : nd, nd + 1 : end)' * B( : ,nd + 1 : end)'));
                        C = [C1; C2];
                end
        end
    end
end
