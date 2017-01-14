function [ W, res0 ] = init_res_dae2_so( eqn, opts, RHS)
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
for mat='MDKG'
    if(~isfield(eqn, sprintf('%c_',mat)) || ~eval(sprintf('isnumeric(eqn.%c_)',mat)))
        error('MESS:error_arguments', 'field eqn.%c_ is not defined',mat);
    end
end

nv = size(eqn.M_,1);
np = size(eqn.G_,1);

if ~isfield(eqn,'type')
  eqn.type='N';
  warning('MESS:equation_type',['Unable to determine type of equation.'...
    'Falling back to type ''N''']);
end

if (~isnumeric(RHS)) || (~ismatrix(RHS))
    error('MESS:error_arguments','RHS has to be a matrix');
end

%% compute low rank residual
RHStemp = zeros(2*nv+np,size(RHS,2));
RHStemp(1:size(RHS,1),:)=RHS;

%if eqn.type=='N'
%    S = [speye(nv,nv),sparse(nv,nv),sparse(nv,np);
%          sparse(nv,nv),eqn.M_,eqn.G_';
%          sparse(np,nv),eqn.G_,sparse(np,np)]; 
%     X = full( S \ RHStemp);
%     W = [X(1:nv,:);eqn.M_*X(nv+1:2*nv,:)];
% else
%     S = [speye(nv,nv),sparse(nv,nv),sparse(nv,np);
%          sparse(nv,nv),eqn.M_',eqn.G_';
%          sparse(np,nv),eqn.G_,sparse(np,np)]; 
%     %X = full( S \ [RHS; sparse(np,size(RHS, 2))]);
%     X = full( S \ RHStemp);
%     W = [X(1:nv,:);eqn.M_'*X(nv+1:2*nv,:)];
% end

S = [speye(nv,nv),sparse(nv,nv),sparse(nv,np);
     sparse(nv,nv),eqn.M_,eqn.G_';
     sparse(np,nv),eqn.G_,sparse(np,np)]; 
X = full( S \ RHStemp);
W = [X(1:nv,:);eqn.M_*X(nv+1:2*nv,:)];



%% compute res0
if opts.adi.LDL_T
    res0 = max(abs(eig(W' * W * eqn.S)));
else
    res0 = norm(W' * W, 2);
end

end

