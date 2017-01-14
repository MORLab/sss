function [eqn,erg] = init_dae_1(eqn, opts, flag1, flag2)
%% init(eqn, flagA, flagE)
% return true or false if Data for A_ and E_ resp. flag1 and flag2  are a
% vailabe and correct in eqn.
%
%   erg = init(eqn,flag1);
%   erg = init(eqn,flag1,flag2);
%
%   erg = init(eqn,'A')    (==init(eqn,'A','A'));
%   erg = init(eqn,'E')    (==init(eqn,'E','E'));
%   erg = init(eqn,'A','E')  (==init(eqn,'E','A'));
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
%   uses no other dae_1 functions

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
    error('MESS:control_data','Number of input Arguments are at least 3');

%% erg = init(eqn, flag1);    
elseif(na==3)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
        case {'E','e'}
            [eqn,erg] = checkE(eqn);
        otherwise
            error('MESS:control_data','flag1 has to be ''A_'' or ''E_''');
    end
    
%% erg = init(eqn,flag1,flag2);
elseif(na==4)
    switch flag1
        case {'A','a'}
            [eqn,erg] = checkA(eqn);
            switch flag2
                case {'A','a'}
                    [eqn,ergA]= checkA(eqn);
                    erg = erg && ergA;
                case {'E','e'}
                    [eqn,ergE]= checkE(eqn);
                    erg = erg && ergE;
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
                    erg =  erg && ergE;
                otherwise
                    error('MESS:control_data','flag2 has to be ''A'' or ''E''');
            end
        otherwise
            error('MESS:control_data','flag1 has to be ''A'' or ''E''');
    end 
end
%% Compute reduced B and C
st = eqn.st;
if size(eqn.B, 1) > st
        eqn.B = eqn.B(1 : st, :) - eqn.A_(1 : st, st + 1 : end) ...
            * (eqn.A_(st + 1 : end, st + 1 : end) \ eqn.B(st + 1 : end, :));
end
if size(eqn.C, 2) > st
        eqn.C = eqn.C( : , 1 : st) - (eqn.C( : , st + 1 : end) ...
            / eqn.A_(st +1 : end, st + 1 : end)') * eqn.A_(1 : st, st+1 : end)';
end
end

%% checkdata for A_
function [eqn,erg] = checkA(eqn)
if ~isfield(eqn,'A_') || ~isnumeric(eqn.A_)
    error('MESS:equation_data',...
      'Empty or Corrupted field A detected in equation structure.')
end
if  (size(eqn.A_,1) ~= size(eqn.A_,2))
    error('MESS:error_arguments', 'field eqn.A_ has to be quadratic');
end
if(~issparse(eqn.A_))
    warning('MESS:control_data','A is not sparse');
end
if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.')
end
erg = 1;
end

%% checkdata for E_
function [eqn,erg] = checkE(eqn)
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.')
elseif (~isfield(eqn,'E_') || ~isnumeric(eqn.E_)) && eqn.haveE
    error('MESS:equation_data',...
      'Empty or Corrupted field E detected in equation structure.')
end
st = eqn.st;
if ~isfield(eqn,'A_') || ~isnumeric(eqn.A_)
    error('MESS:equation_data',...
      'Empty or Corrupted field A detected in equation structure.')
end
n=size(eqn.A_,1);
  
if ~eqn.haveE
  % E = [ I 0 ]
  %     [ 0 0 ]
  eqn.E_=sparse(1:st,1:st,ones(st, 1),n,n,st);
else
  if  (size(eqn.E_,1) ~= size(eqn.E_,2))
    error('MESS:error_arguments', 'field eqn.E_ has to be quadratic');
  end
  if(~issparse(eqn.E_))
    warning('MESS:control_data','E is not sparse');
  end
  % check size(A) == size(E)?
  if (n~=size(eqn.E_,1))
    error('MESS:error_arguments','dimensions of E and A must coincide');
  end
  % E = [ E1 0 ]
  %     [ 0  0 ]
  if full(any([any(eqn.E_(1:st, st + 1:end)), any(eqn.E_(st+1:end,:))]))
    warning('MESS:control_data','E has to be non-zero only in st x st block');
  end
  erg = 1;
  % erg: bool; without 'full()' erg: 1x1 sparse
end
end
