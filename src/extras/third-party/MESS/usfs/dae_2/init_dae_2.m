function [eqn,erg] = init_dae_2(eqn, opts, flag1, flag2)
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
% A = [ A1 -G';
%       G   0]
%
if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.')
end
if ~isfield(eqn,'A_') || ~isnumeric(eqn.A_)
    error('MESS:equation_data',...
      'Empty or Corrupted field A detected in equation structure.')
end
if  (size(eqn.A_,1) ~= size(eqn.A_,2))
    error('MESS:error_arguments', 'field eqn.A_ has to be quadratic');
end
if(~issparse(eqn.A_))
    warning('MESS:check_data','A is not sparse');
end
% check for A12=-A21' or A12=A21' 
asymmG=~any(any(eqn.A_(eqn.st+1:end,1:eqn.st)+eqn.A_(1:eqn.st,eqn.st+1:end)'));
symmG=~any(any(eqn.A_(eqn.st+1:end,1:eqn.st)-eqn.A_(1:eqn.st,eqn.st+1:end)'));
if ~full( asymmG || symmG ) 
   error('MESS:equation_data',...
      'Corrupted field A detected in equation structure. ')   
   
end
% check if lower right block is empty
if (any(any(eqn.A_(eqn.st+1:end,eqn.st+1:end)))) 
   error('MESS:equation_data',...
      'Corrupted field A detected in equation structure.')   
end
erg = 1;
end

%% checkdata for E_
function [eqn,erg] = checkE(eqn)
if ~isfield(eqn, 'haveE'), eqn.haveE = 0; end
if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
          ['Missing or Corrupted st field detected in equation ' ...
           'structure.'])
end
if eqn.haveE
  if ~isfield(eqn,'E_') || ~isnumeric(eqn.E_)
    error('MESS:equation_data',...
      'Empty or Corrupted field E detected in equation structure.')
  end
  if  (size(eqn.E_,1) ~= size(eqn.E_,2))
    error('MESS:error_arguments', 'field eqn.E_ has to be quadratic');
  end
  if(~issparse(eqn.E_))
    warning('MESS:check_data','E is not sparse');
  end
  st = eqn.st;
  % E = [ E1 0;
  %       0  0]
  if full(any([any(eqn.E_(1:st, st + 1:end)), any(eqn.E_(st+1:end,:))]))
    warning('MESS:check_data',['E has to be non-zero only in the ' ...
                        'upper left st x st block']);
  end
else
  % E = [ I 0 ]
  %     [ 0 0 ]
  if ~isfield(eqn,'A_') || ~isnumeric(eqn.A_)
    error('MESS:equation_data',...
      'Empty or Corrupted field A detected in equation structure.')
  end
  st = eqn.st;
  n=size(eqn.A_,1);
  eqn.E_=sparse(1:st,1:st,ones(st, 1),n,n,st);
end
erg = 1;
% erg: bool; without 'full()' erg: 1x1 sparse
end
