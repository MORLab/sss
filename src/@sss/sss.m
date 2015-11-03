classdef sss
% SSS - Sparse state-space LTI system (sss) class
%
% Syntax:
%       sys = SSS(A,B,C)
%       sys = SSS(A,B,C,D,E)
%       sys = SSS(A,B,C,D,E,Ts)
%       sys = SSS(sys_ss)
%       sys = SSS(D)
%
% Description:
%       TODO
%
% Input Arguments:
%       -A: system matrix
%       -B: input matrix
%       -C: output matrix
%       -D: static gain matrix
%       -E: descriptor matrix
%       -Ts: sampling time
%       -sys_ss: control system toolbox state-space (ss)-object
%
% Output Arguments:
%       -sys: sparse state-space (sss)-object
%
% Examples:
%       To create a sparse state-space model of the benchmark "build"
%       (SISO) use:
%
%> load build.mat
%> sys = sss(A,B,C)
%> sysd = diag(sys)
%> [Ad,Bd,Cd] = ssdata(sysd);
%
%       SSS also supports descriptor and DAE systems. To create a sparse
%       state-space model of the benchmark "PEEC_MTLn1600" (DSSS, MIMO)
%       use:
% 
%> load PEEC_MTLn1600.mat
%> sys = sss(A,B,C,zeros(14,14),E)
%
% See Also: 
%        ss, dss
%
% References:
%		* *[1] Documentation of the Control System Toolbox from MATLAB
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order 
% Reduction and System Analysis Toolbox developed at the Chair of 
% Automatic Control, Technische Universitaet Muenchen. For updates 
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Sylvia Cremer, Thomas Emmert (emmert@tfd.mw.tum.de)
%               Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  03 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.
    
    properties
        A,B,C,D,E,x0,Ts
        InputName, StateName, OutputName
        InputGroup, StateGroup, OutputGroup
        InputDelay, InternalDelay, OutputDelay
        UserData
        Name
    end
    properties(Dependent)
    end
    properties(Dependent, Hidden)
        a,b,c,d,e
        n,p,m
        isSiso, isSimo, isMiso, isMimo, isBig
        isDae, isDescriptor
        u,y
    end
    properties(Hidden)
        poles, residues, invariantZeros
        
        H_inf_norm = []
        H_inf_peakfreq = []
        H_2_norm = []
        
        HankelSingularValues
        T_bal, T_bal_inv
        ConGram, ConGramChol
        ObsGram, ObsGramChol
        
        simulationTime
        decayTime
        
        morInfo
        % morInfo must be a struct containing the fields 'time', 'method', 'orgsys'
    end
    
    methods
        function sys = sss(varargin)
            sys.Ts = 0;
            sys.A = [];sys.B = [];sys.C = [];sys.D = [];sys.E = [];sys.x0 = [];
            sys.InputGroup = struct(); sys.StateGroup = struct(); sys.OutputGroup = struct();
            if nargin==0
                % empty state space
                return
            elseif nargin==1
                %% Identity
                if isa (varargin{1}, 'sss')
                    sys = varargin{1};
                    return
                end
                %% Convert matlab LTI systems
                if (isa(varargin{1},'idpoly'))||(isa(varargin{1},'idtf'))
                    varargin{1} = ss(varargin{1},'augmented');
                end
                if isa(varargin{1},'tf')
                    varargin{1} = ss(varargin{1});
                end
                if isa(varargin{1}, 'ss') % convert ss to sss
                    sys_ss = varargin{1};
                    
                    sys.A = sys_ss.A;
                    sys.B = sys_ss.B;
                    sys.C = sys_ss.C;
                    sys.D = sys_ss.D;
                    sys.E = sys_ss.E;
                    sys.Ts = sys_ss.Ts;
                    
                    sys.y = sys_ss.y;
                    sys.u = sys_ss.u;
                    sys.OutputGroup = sys_ss.OutputGroup;
                    sys.InputGroup = sys_ss.InputGroup;
                    sys.UserData= sys_ss.UserData;
                    sys.Name = sys_ss.Name;
                    sys.StateName = sys_ss.StateName;
                    return
                end
                %% Static gain
                if isnumeric(varargin{1})
                    sys.D = varargin{1};
                    return
                end
            elseif nargin==2
                error('Please specify matrices A, B, C.')
            elseif nargin>=3
                sys.A = varargin{1};
                sys.B = varargin{2};
                sys.C = varargin{3};
            end
            
            if nargin>=4
                sys.D = varargin{4};
            end
            
            if nargin>=5
                sys.E = varargin{5};
            end
            
            if nargin >= 6
                sys.Ts = varargin{6};
            end
            
            if nargin > 6
                error('Too many input arguments.')
            end
        end
        
        %% Get Basic Properties
        function m = get.m(sys) % number of inputs
            m = size(sys.B,2);
        end
        function n = get.n(sys) % system order
            n = size(sys.A,1);
        end
        function p = get.p(sys) % number of outputs
            p = size(sys.C,1);
        end
        
        function D = get.D(sys)
            D = sys.D;
            if isempty(D)
                D = sparse(sys.p,sys.m);
            end
        end
        function E = get.E(sys)
            E = sys.E;
            if isempty(E)
                E = speye(sys.n);
            end
        end
        
        function x0 = get.x0(sys)
            x0 = sys.x0;
            if isempty(x0)
                x0 = sparse(1:sys.n,1,0);
            end
        end
        
        %% Set basic properties
        function sys = set.A(sys, A)
            if size(A,1) ~= size(A,2)
                error('A must be square.')
            end
            sys.A = sparse(A);
            sys.poles=[];
        end
        
        function sys = set.B(sys, B)
            if size(B,1) ~= size(sys.A,1)
                error('A and B must have the same number of rows.')
            end
            sys.B = sparse(B);
        end
        
        function sys = set.C(sys, C)
            if size(C,2) ~= size(sys.A,2)
                error('A and C must have the same number of columns.')
            end
            sys.C = sparse(C);
        end
        
        function sys = set.D(sys, D)
            D = sparse(D);
            
            if isempty(D)
                sys.D = [];
            elseif isempty(sys.A)&&isempty(sys.B)&&isempty(sys.C)
                m_ = size(D,2);
                n_ = 0;
                p_ = size(D,1);
                
                sys.A = sparse(n_,n_);
                sys.B = sparse(n_,m_);
                sys.C = sparse(p_,n_);
                sys.D = D;
            else
                if size(D,2) ~= sys.m
                    error('B and D must have the same number of columns.')
                end
                if size(D,1) ~= sys.p
                    error('C and D must have the same number of rows.')
                end
                sys.D = D;
            end
        end
        
        function sys = set.E(sys, E)
            E = sparse(E);
            
            if isempty(E)
                sys.E = [];
            else
                if (size(E) ~= size(sys.A))
                    error('E and A must have the same size.')
                end
                % check whether descriptor matrix is not unity
                if any(any(E-speye(size(E))))
                    sys.E = E;
                else
                    sys.E = [];
                end
            end
            sys.poles=[];
        end
        
        function sys = set.x0(sys, x0)
            if (~isempty(x0)) && (any(size(x0) ~= [sys.n,1]))
                error('A and x0 must have the same number of rows.')
            end
            sys.x0 = sparse(x0);
        end
        
        function sys = set.Ts(sys, Ts)
            if isscalar(Ts)
                sys.Ts = Ts;
            else
                error('Ts must be scalar.')
            end
        end
        
        %% Compatibility with small letters
        function a = get.a(sys); a = sys.A; end
        function b = get.b(sys); b = sys.B; end
        function c = get.c(sys); c = sys.C; end
        function d = get.d(sys); d = sys.D; end
        function e = get.e(sys); e = sys.E; end
        
        function sys = set.a(sys, a); sys.A = a; end
        function sys = set.b(sys, b); sys.B = b; end
        function sys = set.c(sys, c); sys.C = c; end
        function sys = set.d(sys, d); sys.D = d; end
        function sys = set.e(sys, e); sys.E = e; end
        
        %% Get helper functions
        function isSiso = get.isSiso(sys); isSiso=(sys.p==1)&&(sys.m==1); end
        function isSimo = get.isSimo(sys); isSimo=(sys.p>1)&&(sys.m==1); end
        function isMiso = get.isMiso(sys); isMiso=(sys.p==1)&&(sys.m>1); end
        function isMimo = get.isMimo(sys); isMimo=(sys.p>1)||(sys.m>1); end
        
        function isBig = get.isBig(sys); isBig=(sys.n>5000);end
        
        function isDae = get.isDae(sys)
            if condest(sys.E)==Inf
                isDae = 1;
            else
                isDae = 0;
            end
        end
        
        function isDescriptor = get.isDescriptor(sys)
            isDescriptor = logical(full(any(any(sys.E-speye(size(sys.E))))));
        end
        
        function sys = resolveDescriptor(sys)
            sys.A = sys.E\sys.A;
            sys.B = sys.E\sys.B;
            sys.E = [];
        end
        
        % Return system matrices
        function [A,B,C,D,E,Ts] = dssdata(sys)
            A=sys.A; B=sys.B; C=sys.C; D=sys.D; E=sys.E; Ts=sys.Ts;
        end
        function [A,B,C,D,Ts] = ssdata(sys); [A,B,C,D,~,Ts]=dssdata(sys); end
        
        %% Overload Brackets sys.([],[]) to select I/O channels
        function [varargout] = subsref(sys, arg)
            % Returns selected I/O-channel of a sparse LTI MIMO system
            if strcmp(arg(1).type, '()')
                if length(arg(1).subs)==2
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
        
        %% Delay set and get functions
        function sys = set.InputDelay(sys, del); sys.InputDelay = del; end
        function sys = set.OutputDelay(sys, del); sys.OutputDelay = del; end
        function sys = set.InternalDelay(sys, del); sys.InternalDelay = del; end
        
        function del = get.InputDelay(sys); del = sys.InputDelay; end
        function del = get.OutputDelay(sys); del = sys.OutputDelay; end
        function del = get.InternalDelay(sys); del = sys.InternalDelay; end
        
        %% Input, state and output name set and get functions
        % Output
        function yname = get.y(sys); yname = sys.OutputName; end
        function sys = set.y(sys,name); sys.OutputName = name; end
        function name = get.OutputName(sys)
            name = cell(repmat({''}, sys.p, 1));
            name(1:size(sys.OutputName,1),1) = sys.OutputName;
        end
        function sys = set.OutputName(sys, name)
            if (length(name)==sys.p)
                sys.OutputName = name;
            else
                error('Output label vector too long.')
            end
        end
        % State
        function name = get.StateName(sys)
            name = cell(repmat({''}, sys.n, 1));
            name(1:size(sys.StateName,1),1) = sys.StateName;
        end
        function sys = set.StateName(sys, name)
            if (length(name)==sys.n)
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
            name(1:size(sys.InputName,1),1) = sys.InputName;
        end
        function sys = set.InputName(sys, name)
            if (length(name)==sys.m)
                sys.InputName = name;
            else
                error('Input label vector too long.')
            end
        end
    end
    
    methods(Access = private)
        sys = connect_sss(sys, K)
    end
    
    methods(Static)
        [y,x_,index] = simBackwardEuler(A,B,C,D,E,u,x,Ts,TsSample)
        [y,x_,index] = simDiscrete(A,B,C,D,E,u,x,Ts,TsSample,isDescriptor)
        [y,x_,index] = simForwardEuler(A,B,C,D,E,u,x,Ts,TsSample,isDescriptor)
        [y,x_,index] = simRK4(A,B,C,D,E,u,x,Ts,TsSample,isDescriptor)
    end
    
end

