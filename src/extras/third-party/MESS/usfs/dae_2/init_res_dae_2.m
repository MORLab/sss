function [ W, res0 ] = init_res_dae_2( eqn, opts, RHS)
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
%   uses no other dae_2 function

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
if ~isfield(eqn,'A_') || ~isnumeric(eqn.A_)
    error('MESS:equation_data',...
      'Empty or Corrupted field A detected in equation structure.')
end
if ~isfield(eqn, 'st')    || ~isnumeric(eqn.st)
    error('MESS:st',...
    'Missing or Corrupted st field detected in equation structure.')
end
if ~isfield(eqn,'type')
  eqn.type='N';
  warning('MESS:equation_type',['Unable to determine type of equation.'...
    'Falling back to type ''N''']);
end
if ~isfield(eqn,'E_') || ~isnumeric(eqn.E_)
    error('MESS:equation_data',...
      'Empty or Corrupted field E detected in equation structure.')
end
n=size(eqn.A_,1);
if (~isnumeric(RHS)) || (~ismatrix(RHS))
    error('MESS:error_arguments','RHS has to ba a matrix');
end
if (eqn.st ~= size(RHS, 1))
    error('MESS:error_arguments','eqn.st differs with number of rows of RHS');
end

%% compute low rank residual
    

if eqn.type=='N'
    S = eqn.A_;
    S(1:eqn.st,1:eqn.st) = eqn.E_(1:eqn.st,1:eqn.st);
    % S = [ E -J']
    %     [ J  0 ]
    W = eqn.E_ * full( S \ [RHS; sparse(n-eqn.st,size(RHS, 2))]);
    W = W(1:eqn.st,:);
else
    S = eqn.A_';
    S(1:eqn.st,1:eqn.st) = eqn.E_(1:eqn.st,1:eqn.st)';
    % S = [  E' J']
    %     [ -J  0 ]
    W = eqn.E_' * full( S \ [RHS; sparse(n-eqn.st,size(RHS, 2))]);
    W = W(1:eqn.st,:);
end
%% compute res0
if opts.adi.LDL_T
    res0 = max(abs(eig(W' * W * eqn.S)));
else
    res0 = norm(W' * W, 2);
end

end

