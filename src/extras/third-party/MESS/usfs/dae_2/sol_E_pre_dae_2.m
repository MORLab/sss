function [ eqn, opts, oper ] = sol_E_pre_dae_2( eqn, opts, oper )
%% function pre initializes data and/or functions
%
% Input:
%    eqn    struct contains data for equations
%    
%    opts   struct contains parameters for the algorithm
%   
%    oper   struct contains function handles for operation with A
%
% Output:
% eqn
% opts
% oper

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
  alpha=-1/50;
  if isfield(eqn,'st')&&isnumeric(eqn.st)
    st=eqn.st;
  else
    error('MESS:wrong_arguments','missing or corrupted field st detected');
  end
  if ~isfield(eqn,'S_')
    if(~isfield(eqn,'E_') || ~isnumeric(eqn.E_)...
            || ~isfield(eqn,'A_')) || ~isnumeric(eqn.A_)
        error('MESS:error_arguments','field eqn.E_ or eqn.A_ is not defined or corrupted');
    end
    eqn.S_=alpha*eqn.A_;
    eqn.S_(1:st,1:st)=eqn.E_(1:st,1:st);
    eqn.Scount=1;
  else
    if(~isfield(eqn, 'Scount')) || ~isnumeric(eqn.Scount)
        error('MESS:error_arguments', ['field eqn.Scount is not defined. Did ' ...
                        'you forget to run sol_E_pre?']);
    end
    eqn.Scount=eqn.Scount+1;
  end
end

