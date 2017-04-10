classdef sso
% SSO - Create sparse second-order (sso) model, convert to sso model
%
% Syntax:
%       sys = SSO(M,D,K)
%       sys = SSO(M,D,K,B)
%       sys = SSO(M,D,K,B,Cp)
%       sys = SSO(M,D,K,B,Cp,Cv)
%       sys = SSO(M,D,K,B,Cp,Cv,Df)
%       sys = SSO(Df)
%
% Description:
%       This class allows you to create sparse second-order (sso) model objects by
%       just passing the corresponding sparse system matrices.
%
%
% Input Arguments:
%
% Output Arguments:
%
% Examples:
%
% See Also: 
%        sss, ss, dss
%
% References:
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
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  01 Apr 2017
% Copyright (c) 2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    properties
        M,D,K,B,Cp,Cv,Df,x0
        InputName, StateName, OutputName
        InputGroup, StateGroup, OutputGroup
        UserData
        Name
    end
    properties(Dependent)
        isSiso, isSimo, isMiso, isMimo, isBig
        isDae, isDescriptor
        isSym
    end
    properties(Dependent, Hidden)
        n,p,m
        u,y
    end
    properties(Hidden)
                      
    end
    
    methods
        function sys = sso(varargin)
            sys.M = []; sys.D = []; sys.K = [];
            sys.B = [];
            sys.Cp = [];sys.Cv = [];
            sys.Df = []; 
            sys.x0 = [];
            sys.InputGroup = struct(); sys.StateGroup = struct(); sys.OutputGroup = struct();
            if nargin==0
                % empty state space
                return
            elseif nargin==1
                %% Identity
                if isa (varargin{1}, 'sso')
                    sys = varargin{1};
                    return
                end
                %% Static gain
                if isnumeric(varargin{1})
                    sys.Df = varargin{1};
                    return
                end
            elseif nargin==2
                error('Please specify matrices M, D, K.')
            elseif nargin>=3
                sys.M = varargin{1};
                sys.D = varargin{2};
                sys.K = varargin{3};
            end
            
            if nargin>=4
                sys.B = varargin{4};
            end
            
            if nargin>=5
                sys.Cp = varargin{5};
            end
            
            if nargin >= 6
                sys.Cv = varargin{6};
            end
            
            if nargin >= 7
                sys.Df = varargin{7};
            end
            
            if nargin > 7
                error('Too many input arguments.')
            end
        end
        
        %% Get Basic Properties
        function m = get.m(sys) % number of inputs
            m = size(sys.B,2);
        end
        function n = get.n(sys) % system order
            n = size(sys.M,1);
        end
        function p = get.p(sys) % number of outputs
            p = size(sys.Cp,1);
        end
        
        function Df = get.Df(sys)
            Df = sys.Df;
            if isempty(Df)
                Df = sparse(sys.p,sys.m);
            end
        end
        function Cv = get.Cv(sys)
            Cv = sys.Cv;
            if isempty(Cv)
                Cv = sparse(sys.p,sys.n);
            end
        end
        
        function x0 = get.x0(sys)
            x0 = sys.x0;
            if isempty(x0)
                x0 = sparse(1:sys.n,1,0);
            end
        end
        
        %% Set basic properties
        function sys = set.M(sys, M)
            if size(M,1) ~= size(M,2)
                error('M must be square.')
            end
            sys.M = sparse(M);
        end
        
        function sys = set.D(sys, D)
            if size(D) ~= size(sys.M)
                error('M, D, K must have the same size')
            end
            sys.D = sparse(D);
        end
        
        function sys = set.K(sys, K)
            if size(K) ~= size(sys.M)
                error('M, D, K must have the same size')
            end
            sys.K = sparse(K);
        end
        
        function sys = set.B(sys, B)
            if size(B,1) ~= size(sys.M,1)
                error('B and M must have the same number of rows')
            end
            sys.B = sparse(B);
        end
        
        function sys = set.Cp(sys, Cp)
            if size(Cp,2) ~= size(sys.M,2)
                error('M and Cp must have the same number of columns.')
            end
            sys.Cp = sparse(Cp);
        end
        
        function sys = set.Cv(sys, Cv)
            if size(Cv,2) ~= size(sys.M,2)
                error('M and Cv must have the same number of columns.')
            end
            sys.Cv = sparse(Cv);
        end
        
        function sys = set.Df(sys, Df)
            Df = sparse(Df);
            
            if isempty(Df)
                sys.Df = [];
            elseif isempty(sys.M)&&isempty(sys.B)&&isempty(sys.Cp)
                m_ = size(Df,2);
                n_ = 0;
                p_ = size(Df,1);
                
                sys.M   = sparse(n_,n_);
                sys.D   = sparse(n_,n_);
                sys.K   = sparse(n_,n_);
                sys.B   = sparse(n_,m_);
                sys.Cf  = sparse(p_,n_);
                sys.Cv  = sparse(p_,n_);
                sys.Df  = Df;
            else
                if size(Df,2) ~= sys.m
                    error('B and Df must have the same number of columns.')
                end
                if size(Df,1) ~= sys.p
                    error('Cp and Df must have the same number of rows.')
                end
                sys.Df = Df;
            end
        end
        
        
        function sys = set.x0(sys, x0)
            if (~isempty(x0)) && (any(size(x0) ~= [sys.n,1]))
                error('M and x0 must have the same number of rows.')
            end
            sys.x0 = sparse(x0);
        end
        
        function sys = set.isSym(sys, isSym)
            sys.isSym = isSym;
        end
        
        
        %% Get helper functions
        function isSiso = get.isSiso(sys); isSiso=(sys.p==1)&&(sys.m==1); end
        function isSimo = get.isSimo(sys); isSimo=(sys.p>1)&&(sys.m==1); end
        function isMiso = get.isMiso(sys); isMiso=(sys.p==1)&&(sys.m>1); end
        function isMimo = get.isMimo(sys); isMimo=(sys.p>1)||(sys.m>1); end
        
        function isBig = get.isBig(sys); isBig=(sys.n>5000);end
        
        function isDae = get.isDae(sys)
            if condest(sys.M)==Inf
                isDae = true;
            else
                isDae = false;
            end
        end
        
        function isDescriptor = get.isDescriptor(sys)
            isDescriptor = logical(full(any(any(sys.M-speye(size(sys.M))))));
        end
        
        function sys = resolveDescriptor(sys)
            temp = sys.M\[sys.D,sys.K,sys.B]; %only 1 LU decomposition
            sys.D = temp(:,1:sys.n);
            sys.K = temp(:,sys.n+1:2*sys.n);
            sys.B = temp(:,2*sys.n+1:end);
            sys.M = [];
        end
        
        function isSym = get.isSym(sys)
           if isequal(sys.isSym,0) || isequal(sys.isSym,1)
                isSym = sys.isSym;
            else
                if full(max(max(sys.M-sys.M.')))<1e-6 && full(max(max(sys.D-sys.D.')))<1e-6 && full(max(max(sys.K-sys.K.')))<1e-6 
                    isSym = true;
                else
                    isSym = false;
                end
            end
        end
        
        % Return system matrices
        function [M,D,K,B,Cp,Cv,Df] = dssdata(sys)
            M=sys.M; D=sys.D; K=sys.K; B=sys.B; Cp=sys.Cp; Cv=sys.Cv; Df=sys.Df; 
        end
        function [M,D,K,B,Cp,Cv,Df] = ssdata(sys); [M,D,K,B,Cp,Cv,Df]=dssdata(sys); end
        
        % Detect empty sss-models
        function empty = isempty(sys)
            if sys.n == 0
                empty = true;
            else
                empty = false;
            end
        end
        
        %% Overload Brackets sys.([],[]) to select I/O channels
        function [varargout] = subsref(sys, arg)
            % Returns selected I/O-channel of a sparse LTI MIMO system
            if strcmp(arg(1).type, '()')
                if length(arg(1).subs)==2
                    %change ':' to actual indices
                    if strcmp(arg(1).subs{1},':')
                        arg(1).subs{1} = 1:sys.p;
                    end
                    if strcmp(arg(1).subs{2},':')
                        arg(1).subs{2} = 1:sys.m;
                    end
                            
                    sys = sys.truncate(arg(1).subs{1}, arg(1).subs{2});
                    if length(arg)==1
                        varargout = {sys};
                        return
                    end
                end
            elseif strcmp(arg(1).type, '.')
                [varargout{1:nargout}] = builtin('subsref',sys,arg);
                return
            end
            if length(arg)>1
                [varargout{1:nargout}] = builtin('subsref',sys,arg(2:end));
            end
        end
  
        %% Input, state and output name set and get functions
        % Output
        function yname = get.y(sys); yname = sys.OutputName; end
        function sys = set.y(sys,name); sys.OutputName = name; end
        function name = get.OutputName(sys)
            name = cell(repmat({''}, sys.p, 1));
            if ~isempty(sys.OutputName)
                name(1:size(sys.OutputName,1),1) = sys.OutputName;
            end
        end
        function sys = set.OutputName(sys, name)
            if isempty(name) || all(cellfun(@isempty, name))
                sys.OutputName = [];
            elseif (length(name)==sys.p)
                sys.OutputName = name;
            else
                error('Output label vector too long.')
            end
        end
        % State
        function name = get.StateName(sys)
            name = cell(repmat({''}, sys.n, 1));
            if ~isempty(sys.StateName)
                name(1:size(sys.StateName,1),1) = sys.StateName;
            end
        end
        function sys = set.StateName(sys, name)
            if isempty(name) || all(cellfun(@isempty, name))
                sys.StateName = [];
            elseif (length(name)==sys.n)
                sys.StateName = name;
            else
                error('State label vector too long.')
            end
        end
        % Input
        function name = get.u(sys); name = sys.InputName; end
        function sys = set.u(sys,name); sys.InputName = name; end
        function name = get.InputName(sys)
            name = cell(repmat({''}, sys.m, 1));
            if ~isempty(sys.InputName)
                name(1:size(sys.InputName,1),1) = sys.InputName;
            end
        end
        function sys = set.InputName(sys, name)
            if  isempty(name) || all(cellfun(@isempty, name))
                sys.InputName = [];
            elseif (length(name)==sys.m)
                sys.InputName = name;
            else
                error('Input label vector too long.')
            end
        end
    end
      
    
end

