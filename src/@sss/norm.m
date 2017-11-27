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
%       The $H_\infty$ norm of a system is computed using a Newton Method
%       optimization. Firstly, many frequency responses are computed, then, the
%       maximal found frequency response is locally optimized in order to find the $H_\infty$ norm.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: sss-object containing the LTI system
%       *Optional Input Arguments:* 
%       -p: choice of H_2-norm or H_infty-norm 
%           [{'2'} / 'inf']
%       -Opts:              a structure containing following options
%           -.lyapchol:     choose lyapunov equation solver
%                           [{'auto'} / 'adi' / 'hammarling' / 'crksm']
%           -.lse:          solve linear system of equations
%                           [{'sparse'} / 'full' / 'gauss' / 'hess' / 'iterative']
%           -.stabcheck:    perform a stability check
%                           [{true},false]
%
% Output Arguments:
%       -nrm:             value of norm
%       -hInfPeakfreq:    peak frequency of magnitude of H_\infty norm
%
% Examples:
%       The following code computes the $H_2$- and the $H_\infty$-norm of the
%       benchmark 'CDplayer' (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> h2Norm=norm(sys,2)
%> h_infNorm=norm(sys,inf)
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
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Jorge Luiz Moreira Silva, Heiko Panzer, Sylvia Cremer, Rudy Eid
%               Alessandro Castagnotto, Maria Cruz Varona, Lisa Jeschek
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  15 Apr 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------
%%  Define execution parameters
Def.lyapchol    = 'auto'; 
Def.lse         = 'sparse'; %lse 
Def.stabcheck   = true;

%% Computation
if isempty(sys)
    nrm = 0;
    if nargout>1
        varargout{1}=0; %hInfPeakfreq
    end
else
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
        
        if Opts.stabcheck && isstable(sys)~=1
            nrm=Inf;
            return;
        end

        Opts.method=Opts.lyapchol; %translate option for lyapchol function
        try
            R=lyapchol(sys,Opts);
            nrm=norm(sys.C*R,'fro');
        catch ex
            warning(ex.identifier, 'Error solving Lyapunov equation. Trying without Cholesky factorization...')
            X = lyap(sys.A, sys.B*sys.B', [], sys.E);
            nrm=sqrt(trace(sys.C*X*sys.C'));
        end

        if imag(nrm)~=0
            nrm=Inf;
        end
    else
        error(['H_' num2str(p) '-norm not implemented.'])
    end
end
end

function [ H_Infty,freq ] = H_Infty( sys )
 
    if sys.n < 500
        [H_Infty,freq] = norm(ss(sys),inf); %built in
    else

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
        w=freq;
        [Deriv0,Deriv1,Deriv2]=computeDerivatives(A,B,C,D,E,w);
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
            [Deriv0,Deriv1,Deriv2]=computeDerivatives(A,B,C,D,E,w);

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
end

function [Deriv0,Deriv1,Deriv2]=computeDerivatives(A,B,C,D,E,w)
Opts.krylov='standardKrylov';

% tangential directions
Rt=zeros(size(B,2),size(B,2)*3);
for nB2=1:size(B,2)
    Rt(nB2,(nB2-1)*3+1:(nB2-1)*3+3)=ones(1,3);
end

linSolve=solveLse(A,B,E,-ones(1,3*size(B,2))*w*1i, Rt, Opts);

% sort solutions
linSolve1=zeros(size(B));
linSolve2=zeros(size(B));
linSolve3=zeros(size(B));

for nB2=1:size(B,2)
    linSolve1(:,nB2)=linSolve(:,(nB2-1)*3+1);
    linSolve2(:,nB2)=linSolve(:,(nB2-1)*3+2);
    linSolve3(:,nB2)=linSolve(:,(nB2-1)*3+3);
end

resp=(C*linSolve1)+D;
Deriv0=((resp)'*resp); %it is ln
%Compute first derivative
respp=-1i*C*linSolve2;
Deriv1=(respp'*resp)+(respp'*resp)';
%Compute second derivative
resppp=-2*C*linSolve3;
Deriv2=(resppp'*resp)+(resppp'*resp)'+2*(respp'*respp);
end

