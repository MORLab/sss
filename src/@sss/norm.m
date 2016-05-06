function [nrm, varargout] = norm(sys, varargin)
% NORM - Computes the p-norm of an sss LTI system
%
% Syntax:
%       nrm = NORM(sys)
%       nrm = NORM(sys,p)
%       [nrm, hInfPeakfreq] = NORM(sys, inf)
%       nrm = NORM(...,Opts)
%
% Description:
%       This function computes the p-norm of an LTI system given
%       as a sparse state-space (sss) object sys. The value of p can be
%       passed as a second optional argument to the function and is set to
%       2 otherwise.
%       The H_Infinity norm of a system is computed using a Newton Method
%       optimization. Firstly, many frequency responses are computed, then, the
%       maximal found frequency response is locally optimized in order to find the H_Infinity norm.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: sss-object containing the LTI system
%       *Optional Input Arguments:* 
%       -p: choice of H_2-norm or H_inf-norm 
%           [{'2'} / 'inf']
%       -Opts:              a structure containing following options
%           -.adi:          try only solution by adi or lyapunov equation
%                           [{'0'} / 'adi' / 'lyap']
%           -.lse:          solve linear system of equations (only for adi)
%                           [{'gauss'} / 'luChol']
%
% Output Arguments:
%       -nrm:             value of norm
%       -hInfPeakfreq:    peak frequency of magnitude of H_inf norm
%
% Examples:
%       The following code computes the H2- and the H_inf-norm of the
%       benchmark 'CDplayer' (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> h2Norm=norm(sys,2);
%> h_infNorm=norm(sys,inf);
%
% See also:
%       norm, sss, lyapchol
%
% References:
%       * *[1] Antoulas (2005)*, Approximation of large-scale Dynamical Systems
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
% Authors:      Jorge Luiz Moreira Silva, Heiko Panzer, Sylvia Cremer, Rudy Eid
%               Alessandro Castagnotto, Maria Cruz Varona, Lisa Jeschek
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  15 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
%%  Define execution parameters
Def.adi= 0; %use only adi or lyapunov equation ('0','adi','lyap')
Def.lse= 'gauss'; %lse (used only for adi)

p=2;    % default: H_2
if nargin>1
    if isa(varargin{1}, 'double')
        p=varargin{1};
    elseif strcmpi(varargin{1},'inf')
        p=inf;
    elseif isa(varargin{1},'struct')
        Opts=varargin{1};
    else
        error('Input must be ''double''.');
    end
    if nargin==3
        Opts=varargin{2};
    end
end

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

if isinf(p)
    % H_inf-norm
    [nrm, hInfPeakfreq] = H_Infty(sys);
    
    if nargout>1
        varargout{1}=hInfPeakfreq;
    end
    
elseif p==2
    % when D ~=0, then H2 norm is unbounded
    if any(any(sys.D))
        nrm=inf;
        return
    end
    
    % see if a Gramian or its Cholesky factor is already available
    if isempty(sys.ConGramChol)
        if isempty(sys.ObsGramChol)
            if isempty(sys.ConGram)
                if isempty(sys.ObsGram)
                    % No, it is not. Solve Lyapunov equation.
                    if ~sys.isDae
                        % options for mess
                        % eqn struct: system data
                        eqn=struct('A_',sys.A,'E_',sys.E,'B',sys.B,'C',sys.C,'type','N','haveE',sys.isDescriptor);

                        % opts struct: mess options
                        messOpts.adi=struct('shifts',struct('l0',20,'kp',50,'km',25,'b0',ones(sys.n,1),...
                            'info',0),'maxiter',300,'restol',0,'rctol',1e-12,...
                            'info',0,'norm','fro');
                        
                        % user functions: default
                        if strcmp(Opts.lse,'gauss')
                            oper = operatormanager('default');
                        elseif strcmp(Opts.lse,'luChol')
                            if sys.isSym
                                oper = operatormanager('chol');
                            else
                                oper = operatormanager('lu');
                            end
                        end
                    end
                    try
                        if strcmp(Opts.adi,'lyap') || sys.n<100 || sys.isDae || (~strcmp(Opts.adi,'adi') && sys.n<500)
                            error('lyap');
                        end
                        if isstable(sys)~=1
                            warning('System appears to be unstable. The norm will be set to Inf.');
                            nrm=Inf;
                            return;
                        end

                        % get adi shifts
                        [messOpts.adi.shifts.p, eqn]=mess_para(eqn,messOpts,oper);

                        % low rank adi
                        [R,~,eqn]=mess_lradi(eqn,messOpts,oper);

                        nrm=norm(R'*eqn.C','fro');
                    catch ex
                        if strcmp(Opts.adi,'adi')
                            if strcmp(ex.message,'lyap')
                                error('ADI failed. System may be too small or DAE.');
                            else
                                error(ex.message);
                            end
                        elseif  ~strcmp(ex.message,'lyap');
                            warning([ex.message,' Trying without ADI...']);
                        end
                        if sys.isDescriptor
                            try
                                try
                                    sys.ConGramChol = lyapchol(sys.A,sys.B,sys.E); % P=S'*S3
                                    nrm=norm(sys.ConGramChol*sys.C','fro');
                                    if ~isreal(nrm)
                                        error('Gramian must be positive definite');
                                    end
                                catch ex3
                                    P = lyapchol(sys.A',sys.C',sys.E');
                                    nrm=norm(P*sys.B,'fro');
                                end
                            catch ex
                                warning(ex.identifier, 'Error solving Lyapunov equation. Trying without Cholesky factorization...')
                                try
                                    try
                                        X = lyap(sys.A, sys.B*sys.B', [], sys.E);
                                        nrm=sqrt(trace(sys.C*X*sys.C'));
                                        if ~isreal(nrm)
                                            error('Gramian must be positive definite');
                                        end
                                    catch ex3
                                        Y = lyap(sys.A', sys.C'*sys.C, [], sys.E');
                                        nrm=sqrt(trace(sys.B'*Y*sys.B));
                                    end
                                catch ex2
                                    warning(ex2.message, 'Error solving Lyapunov equation. Premultiplying by E^(-1)...')
                                    tmp = sys.E\sys.B;
                                    X = lyap(sys.E\sys.A, tmp*tmp');
                                    nrm=sqrt(trace(sys.C*X*sys.C'));
                                end
                            end
                        else
                            try
                                sys.ConGramChol = lyapchol(sys.A,sys.B);
                                nrm=norm(sys.ConGramChol*sys.C','fro');
                            catch ex
                                if strcmp(ex.identifier,'Control:foundation:LyapChol4');
                                    %Unstable system. Set the norm to infinity
                                    warning('System appears to be unstable. The norm will be set to Inf.')
                                    nrm = Inf;
                                else
                                    warning(ex.message, 'Error solving Lyapunov equation. Trying without Cholesky factorization...')
                                    sys.ConGram = lyap(sys.A, sys.B*sys.B');                
                                    nrm=sqrt(trace(sys.C*sys.ConGram*sys.C'));
                                end
                            end
                        end
                    end
                else
                    nrm=sqrt(trace(sys.B'*sys.ObsGram*sys.B));
                end
            else
                nrm=sqrt(trace(sys.C*sys.ConGram*sys.C'));
            end
        else
            nrm=norm(sys.ObsGramChol*sys.B, 'fro');
        end
    else
        nrm=norm(sys.ConGramChol*sys.C','fro');
    end
    
    if imag(nrm)~=0
        nrm=Inf;
    end
else
    error(['H_' num2str(p) '-norm not implemented.'])
end
end

function [ H_Infty,freq ] = H_Infty( sys )

[mag,w] = freqresp( sys ); %Compute the frequency response of the system
[magExtreme]=freqresp(sys,[0,inf]);
mag=cat(3,magExtreme(:,:,1),mag,magExtreme(:,:,2));
w=[0;w;inf];
possibleIndex=1:numel(w);
%%Defining Boundaries
maxNorm=sqrt(sum(sum(abs(mag).^2,1),2)); %Frobenius Norm always greater or equal than the 2-norm
[~,indexMaxNorm]=max(maxNorm);
index=indexMaxNorm;
H_Infty=norm(mag(:,:,index),2);
i=0;
%Find the maximum norm computed in freqresp
while(1==1)
    i=i+1;
    maxNorm(index)=0;
    possibleIndex(maxNorm<=H_Infty)=[]; %All norms that can't be greater than the value of the variable of H_Infty are discarded
    maxNorm(maxNorm<=H_Infty)=[];
    if not(numel(maxNorm))
        break;
    end
    [~,index]=max(maxNorm);
    computedNorm=norm(mag(:,:,possibleIndex(index)),2);
    if (computedNorm>H_Infty)
        H_Infty=computedNorm;
        indexMaxNorm=possibleIndex(index);
    end
end
freq=w(indexMaxNorm);

%Optimization of the maximum Hinfty norm using newton method
[A,B,C,D,E]=dssdata(sys);
minusA=-A;
w=freq;
[Deriv0,Deriv1,Deriv2]=computeDerivatives(minusA,B,C,D,E,w);
[eigenVectors,eigenValues]=(eig(full(Deriv0)));
eigenValues=diag(eigenValues);
[~,Index]=max(eigenValues);
vec=eigenVectors(:,Index);
firstDeriv=real(vec'*Deriv1*vec); %Computation of first derivative
secondDeriv=real(vec'*Deriv2*vec); %Computation of second derivative
delta=inf;
i=0;
while(1==1) %Newton-Method Iteration
    i=i+1;
    deltaBefore=delta;
    delta=firstDeriv/secondDeriv;
    w=w-delta; %Update of w in order to find firstDeriv=0
    [Deriv0,Deriv1,Deriv2]=computeDerivatives(minusA,B,C,D,E,w);
    
    [eigenVectors,eigenValues]=(eig(full(Deriv0)));
    eigenValues=diag(eigenValues);
    [~,Index]=max(eigenValues);
    vec=eigenVectors(:,Index);
    if (abs((norm(deltaBefore)-norm(delta))/norm(delta))<1e-9)|norm(delta(end))<eps(w)|i>=10 %Condition to stop iterations
        break;
    end
    firstDeriv=real(vec'*Deriv1*vec);
    secondDeriv=real(vec'*Deriv2*vec);
end
freq=w;
H_Infty=sqrt(max(eig(full(Deriv0))));
end

function [Deriv0,Deriv1,Deriv2]=computeDerivatives(minusA,B,C,D,E,w)
[L,U,k,l,S]=lu((minusA+E*w*1i),'vector');
b=S\B; b=b(k,:);
LinearSolve0=L\b;
LinearSolve0(l,:)=U\LinearSolve0;
resp=(C*LinearSolve0)+D;
Deriv0=((resp)'*resp); %it is ln
%Compute first derivative
b=S\(E*LinearSolve0); b=b(k,:);
LinearSolve1=L\b;
LinearSolve1(l,:)=U\LinearSolve1;
respp=-1i*C*LinearSolve1;
Deriv1=(respp'*resp)+(respp'*resp)';
%Compute second derivative
b=S\(E*LinearSolve1); b=b(k,:);
LinearSolve2=L\b;
LinearSolve2(l,:)=U\LinearSolve2;
resppp=-2*C*LinearSolve2;
Deriv2=(resppp'*resp)+(resppp'*resp)'+2*(respp'*respp);
end

