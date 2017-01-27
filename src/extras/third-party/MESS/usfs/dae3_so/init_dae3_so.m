function [eqn,erg] = init_dae3_so(eqn, opts, flag1, flag2)
%% init_dae_2(eqn, flagA, flagE)
% return true or false if Data for A_ and E_ resp. flag1 and flag2  are 
% availabe and correct in eqn.
%
%   erg = init_dae_2(eqn,flag1);
%   erg = init_dae_2(eqn,flag1,flag2);
%
%   erg = init_dae_2(eqn,'A')    (==init_dae_2(eqn,'A','A'));
%   erg = init_dae_2(eqn,'E')    (==init_dae_2(eqn,'E','E'));
%   erg = init_dae_2(eqn,'A','E')  (==init_dae_2(eqn,'E','A'));
%
%   Input:
%   flag1/flag2    'A'/'E' for checking A or E in eqn
%   eqn             structure with the data
%   opts            struct contains parameters for the algorithm
%
%   Output:
%   eqn             structure with the data
%   opts            struct contains parameters for the algorithm

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
na = nargin;
if(na<=2)
    error('MESS:check_data','Number of input Arguments must be at least 3');
    
    %% erg = init_dae_2(eqn, flag1);
elseif(na==3)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
        case {'E','e'}
            [eqn,erg] = checkE(eqn);
        otherwise
            error('MESS:check_data','flag1 has to be ''A'' or ''E''');
    end
    
    %% erg = init_dae_2(eqn,flag1,flag2);
elseif(na==4)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn,ergE] = checkE(eqn);
                    erg = erg &&ergE;
                otherwise
                    error('MESS:check_data','flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            [eqn, erg] = checkE(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn,ergE] = checkE(eqn);
                    erg =  erg && ergE;
                otherwise
                    error('MESS:check_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:check_data','flag1 has to be ''A'' or ''E''');
    end
end
end

%% checkdata for A_
function [eqn,erg] = checkA(eqn)

if ~isfield(eqn,'D_') || ~isnumeric(eqn.D_)
    error('MESS:equation_data',...
        'Empty or Corrupted field D detected in equation structure.')
elseif ~issparse(eqn.D_)
    warning('MESS:control_data','D is not sparse');
end

if ~isfield(eqn,'K_') || ~isnumeric(eqn.K_)
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
elseif ~issparse(eqn.K_)
    warning('MESS:control_data','K is not sparse');
end

if ~isfield(eqn,'G_') || ~isnumeric(eqn.G_)
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
elseif ~issparse(eqn.G_)
    warning('MESS:control_data','G is not sparse');
end

if  (size(eqn.D_,1) ~= size(eqn.D_,2))
    error('MESS:error_arguments', 'field eqn.D_ has to be quadratic');
end

if  (size(eqn.K_,1) ~= size(eqn.K_,2))
    error('MESS:error_arguments', 'field eqn.K_ has to be quadratic');
end

if  (size(eqn.D_,1) ~= size(eqn.G_,2))
    error('MESS:error_arguments', 'field eqn.G_ has invalid number of columns');
end

erg = 1;
end

%% checkdata for E_
function [eqn,erg] = checkE(eqn)
if ~isfield(eqn, 'haveE'), eqn.haveE = 1; end

if ~isfield(eqn,'M_') || ~isnumeric(eqn.M_)
    error('MESS:equation_data',...
        'Empty or Corrupted field M detected in equation structure.')
end
if  (size(eqn.M_,1) ~= size(eqn.M_,2))
    error('MESS:error_arguments', 'field eqn.M_ has to be quadratic');
end
if(~issparse(eqn.M_))
    warning('MESS:check_data','M is not sparse');
end

if ~isfield(eqn, 'alpha') || ~isnumeric(eqn.alpha)
   error('MESS:equation_data',...
         'No parameter alpha given for shifting infinite eigenvalues of the pencil'); 
end
erg=1;
end
