function [ RHS, res0 ] = init_res_dae_1( eqn, opts, RHS)
%% function init_res initializes the low rank residual W and res0
% 
%  Input: 
%     eqn     structure contains data for A, B and C
%     
%     opts    structure contains parameters for the algorithm
%
%     RHS     right hand side matrix
%     
%  Output:
%  
%    W 
%    res0
%
%   uses no other dae_1 function

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

%% check data
if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.')
end
if ~isfield(eqn,'type')
  eqn.type='N';
  warning('MESS:equation_type',['Unable to determine type of equation.'...
    'Falling back to type ''N''']);
end
if ~isfield(eqn,'A_') || ~isnumeric(eqn.A_)
    error('MESS:equation_data',...
      'Empty or Corrupted field A detected in equation structure.')
end
if (~isnumeric(RHS)) || (~ismatrix(RHS))
    error('MESS:error_arguments','RHS has to ba a matrix');
end
if (eqn.st ~= size(RHS, 1))
    error('MESS:error_arguments','number of rows of A_ differs with number of rows of RHS');
end
    
%% compute res0
if opts.adi.LDL_T
    res0 = max(abs(eig(RHS' * RHS * eqn.S)));
else
    res0 = norm(RHS' * RHS, 2);
end

end

