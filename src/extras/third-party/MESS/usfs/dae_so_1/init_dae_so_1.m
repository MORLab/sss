function [eqn,erg] = init_dae_so_1(eqn, opts, flag1, flag2)
%% init_dae_so_1(eqn, flagA, flagE)
% return true or false if Data for A and E resp. flag1 and flag2  are availabe
% and correct in eqn.
%
%   erg = init_so_1(eqn,flag1);
%   erg = init_so_1(eqn,flag1,flag2);
%
%   erg = init_so_1(eqn,'A')    (==init_so_1(eqn,'A','A'));
%   erg = init_so_1(eqn,'E')    (==init_so_1(eqn,'E','E'));
%   erg = init_so_1(eqn,'A','E')  (==init_so_1(eqn,'E','A'));
%
%   Input:
%   flag1/flag2    'A'/'E' for checking A or E in eqn
%   eqn             structure with the data
%   opts            struct contains parameters for the algorithm
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

%% start checking
na = nargin;
if(na<=2)
    error('MESS:control_data','Number of input Arguments are at least 3');

%% erg = init_so_1(eqn, flag1);    
elseif(na==3)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
        case {'E','e'}
            [eqn,erg] = checkE(eqn);
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end
    
%% erg = init_so_1(eqn,flag1,flag2);
elseif(na==4)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn, ergE] = checkE(eqn);
                    erg =  erg && ergE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        case {'E','e'}
            [eqn,erg] = checkE(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA] = checkA(eqn);
                    erg =  erg && ergA;
                case {'E','e'}
                    [eqn,ergE] = checkE(eqn);
                    erg = erg && ergE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end 
end
end

%% checkdata for A
function [eqn,erg] = checkA(eqn)
if (~isfield(eqn,'K_') || ~isnumeric(eqn.K_))
    error('MESS:equation_data',...
        'Empty or Corrupted field K detected in equation structure.')
end
if(~issparse(eqn.K_))
    warning('MESS:control_data','K is not sparse');
end
[n1k, n2k] = size(eqn.K_);
if n1k ~= n2k
    error('MESS:equation_data',...
        'K has to be quadratic')
end
erg = 1;

end

%% checkdata for E
function [eqn,erg] = checkE(eqn)
if (~isfield(eqn,'M_') || ~isnumeric(eqn.M_))
    error('MESS:equation_data',...
        'Empty or Corrupted field M detected in equation structure.')
elseif (~isfield(eqn,'D_') || ~isnumeric(eqn.D_))
    error('MESS:equation_data',...
        'Empty or Corrupted field D detected in equation structure.')
end
if ~isfield(eqn, 'nd')    || ~isnumeric(eqn.nd)
    error('MESS:nd',...
    'Missing or Corrupted nd field detected in equation structure.');
end
if(~issparse(eqn.M_))
    warning('MESS:control_data','M is not sparse');
end
if(~issparse(eqn.D_))
    warning('MESS:control_data','D is not sparse');
end
nd = eqn.nd; 
[n1m, n2m] = size(eqn.M_);
[n1d, n2d] = size(eqn.D_);
if n1m ~= n2m
    error('MESS:equation_data',...
        'M has to be quadratic')
end
if n1d ~= n2d
    error('MESS:equation_data',...
        'D has to be quadratic')
end
if n1m ~= n1d
    error('MESS:equation_data',...
        'M and D must have same size')
end
if full(any([any(eqn.M_(1:nd, nd + 1:end)), any(eqn.M_(nd+1:end,:))]))
    warning('MESS:control_data','M has to be non-zero only in nd x nd block');
end
if full(any([any(eqn.D_(1:nd, nd + 1:end)), any(eqn.D_(nd+1:end,:))]))
    warning('MESS:control_data','D has to be non-zero only in nd x nd block');
end
erg = 1;

end
