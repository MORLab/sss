function varargout = second2first(varargin)
%SECOND2FIRST - Convert a 2nd order system to state space representation
%
% Syntax:
%       sys         = SECOND2FIRST(sysSO)
%       sys         = SECOND2FIRST(M,D,K,B,Cx)
%       sys         = SECOND2FIRST(M,D,K,B,[],Cv)
%       sys         = SECOND2FIRST(M,D,K,B,Cx,Cv)
%       sys         = SECOND2FIRST(M,D,K,B,Cx,Cv,Df,Opts)
%       [A,B,C,D,E] = SECOND2FIRST(M,D,K,B,Cx,Cv,Df,Opts)
%
% Description:
%       This function takes the system matrices of a 2nd order system and
%       converts this system to state space representation. 
%       
%       The representation of a 2nd order system thereby is defined as
%       follows:
%
%       $M  ~ \ddot x + D ~ \dot x + K ~ x = B ~ u$
%
%       $y = C_x ~ x + C_v ~ \dot x$
%
%       State space representation uses the following form:
%
%       <<img/second2first.png>>
%
%       The conversion from second order to first order allows some freedom
%       in the choice of the matrix F. Which value from the possible
%       options should be used to create the matrix F can be specified by
%       passing a option-structure with the desired choice as an additional
%       input argument to the function. 
%
%       It is also possible to specify user defined values for the matrix
%       F. The first possibility to do this is to set the field 
%       *Opts.transf2nd* to a scalar value unequal to zero. This value then
%       is multiplied with the identity matrix to compute F.
%
%       $F = \alpha \cdot I$
%
%       The second possibility is to set the field *Opts.transf2nd* to a
%       positive definite matrix of the correct dimension. Notice that a  
%       suitable choise of the matrix F can improve the condition numbers 
%       of A and E.
%
%       If the representation of the system to convert contains some
%       matrices which are all zero (i.e. no damping), then the value [] 
%       can be specified for these matrices.
%
%       If only one ouput is specified, then SECOND2FIRST return an
%       sss-object. Otherwise, it returns the system matrices of the
%       respective first order model.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -M:     Mass matrix of the second order system
%       -D:     Damping matrix of the second order system
%       -K:     Stiffness matrix of the second order system
%       -B:     Input matrix of the second order system
%       -Cx:    Output matrix of the second order system 
%       *Optional Input Arguments:*
%       -Cv:    Output matrix of the second order system
%       -Opts:	a structure containing following options
%			-transf2nd:     Type of transformation from 2nd to
%                           1st order form: 
%                           [{'I'}, 'K', '-K', 'alpha', scalar , matrix]
%
% Output Arguments:
%       -sys:       Sparse State Space Model (1st order form)
%       -A,B,C,D,E: System matrices in 1st order form
%
% Examples:
%       The following code converts a randomly created 2nd order system
%       without damping to a 1st order system:
%> clear all;
%> sys = second2first(rand(20),[],rand(20),rand(20,1),rand(1,20));
%	
% See Also:
%       loadSss
%
% References:
%       * *[1] B. Salimbahrami*: "Structure Preserving Order Reduction of Large
%       Scale Second Order Models", PhD Thesis at TUM, p. 36f
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Rudy Eid, Niklas Kochdumper,
%               Alessandro Castagnotto
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  10 Nov 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Parse inputs
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end

if isa(varargin{1},'sso')
    [M,D,K,B,Cp,Cv,Df] = dssdata(varargin{1});
else
    M = varargin{1};
    D = varargin{2};
    K = varargin{3};
    if length(varargin)>=4
        B = varargin{4};
    end
    if length(varargin)>=5
        Cp = varargin{5};
    end
    if length(varargin)>=6
        Cv = varargin{6};
    end
    if length(varargin)>=7
        Df = varargin{7};
    end
end
%%
% check wheather options are specified or not 
if ~exist('Opts','var')
    Opts = struct(); 
end

Def = struct('transf2nd','I');
Opts = parseOpts(Opts,Def);


% check the matrix-dimensions of M and D
if size(M,1)~=size(M,2)
    error('M must be symmetric.')
elseif any(size(M)-size(K))
    error('M and K must have same size.');
end

if ~exist('D','var') || isempty(D)
    D = sparse(size(M,1), size(M,1));
elseif any(size(D)-size(K))
    error('D must have same size as M and K.');
end

if ~exist('Df','var') || isempty(Df)
    Df = sparse(size(Cp,1), size(B,2));
elseif any(size(Df)-[size(Cp,1),size(B,2)])
    error('Df must have compatible size with B and C');
end

% check if user defined scalar or matrix is valid
if isnumeric(Opts.transf2nd)
   if ~isreal(Opts.transf2nd)
          error('Value for Opts.transf2nd has to be real.')
   end
   if isscalar(Opts.transf2nd)  % scalar
      if Opts.transf2nd == 0
          error('Value 0 is not allowed for Opts.transf2nd');
      end
   else                         % matrix
      if ~isequal(size(Opts.transf2nd),size(M))
          error('The matrix specified for Opts.transf2nd must have dimesion %d x %d', ...
                size(M,1),size(M,2));
      elseif ~ispd(Opts.transf2nd)
          error('The matrix specified for Opts.transf2nd has to be positive definite.');
      end
   end
end

% evaluate the specified option
if isnumeric(Opts.transf2nd)
    if isscalar(Opts.transf2nd)
        F = Opts.transf2nd * speye(size(K));
    else
        F = Opts.transf2nd;
    end 
else
    switch Opts.transf2nd
        case 'I'
            F = speye(size(K));
        case 'K'
            F = K;
        case '-K'
            F = -K;
        case 'alpha'
            % TODO: add the transformation to strictly dissipative form by
            % Panzer
            warning('alpha option not implemented yet, using K instead');
            alpha = 1; 
            F = alpha*K;                
    end
end

% create the matrices A and E for the first-order-system
n = size(M,1);
O = sparse(n,n);

E = [F O; O M];
A = [O F; -K -D];

% check the matrix-dimensions of Cx and Cv and B
if isempty(Cp) && isempty(Cv)
    % both output vectors are empty
    error('All outputs are zero.');
elseif isempty(Cp)
    % only Cv is given
    Cp = sparse(size(Cv,1),n);
elseif ~exist('Cv', 'var') || isempty(Cv)
    % only Cx is given
    Cv = sparse(size(Cp,1),n);
elseif any(size(Cp)-size(Cv))
    error('Cx and Cv must have same size.');
elseif size(Cp,2)~=size(M,2)
    error('Cx, Cv must have same column dimension as M, D, K')
end

if isempty(B)
    error('All inputs zero.');
elseif size(B,1)~=size(M,1)
    error('B must have same row dimension as M, D, K')
end

% create the matrices C, D and B for the first-order-system 
B = [sparse(n,size(B,2)); B];
C = [Cp, Cv];
D = Df;

if nargout == 1
    % create the first-order-system
    varargout{1}    = sss(A, B, C, D, E);
else
    % return system matrices
    varargout       = {A,B,C,D,E};
end
