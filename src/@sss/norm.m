function [nrm, varargout] = norm(sys, varargin)
% NORM - Computes the p-norm of an sss LTI system
%
% Syntax:
%       nrm = NORM(sys)
%       nrm = NORM(sys,p)
%       [nrm, hInfPeakfreq] = NORM(sys, inf)
%
% Description:
%       This function computes the p-norm of an LTI system given 
%       as a sparse state-space (sss) object sys. The value of p can be 
%       passed as a second optional argument to the function and is set to
%       2 otherwise.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: sss-object containing the LTI system
%       *Optional Input Arguments:* 
%       -p: 2 for H_2-norm or inf for H_inf-norm
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
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Jorge Luiz Moreira Silva, Heiko Panzer, Sylvia Cremer, Rudy Eid
%               Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Dez 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

p=2;    % default: H_2
if nargin>1
    if isa(varargin{1}, 'double')
        p=varargin{1};
    elseif strcmpi(varargin{1},'inf')
        p=inf;
    else
        error('Input must be ''double''.');
    end
end

if isinf(p)
    % H_inf-norm
    if isempty(sys.hInfNorm)
        %warning('calling MATLAB''s built-in norm');
        [sys.hInfNorm, sys.hInfPeakfreq] = H_Infty(sys);
        if nargout>1
            varargout{1}=sys.hInfPeakfreq;
        end
        if inputname(1)
            assignin('caller', inputname(1), sys);
        end
    end
    nrm=sys.hInfNorm; 
elseif p==2
    % H_2-norm
    if ~isempty(sys.h2Norm)
        nrm=sys.h2Norm;
        return
    end
    % wenn D ~=0 ist H2 norm unendlich groﬂ
    if any(any(sys.D))
        nrm=inf;
        sys.h2Norm=inf;
        return
    end

    % see if a Gramian or its Cholesky factor is already available
    if isempty(sys.ConGramChol)
        if isempty(sys.ObsGramChol)
            if isempty(sys.ConGram)
                if isempty(sys.ObsGram)
                    % No, it is not. Solve Lyapunov equation.
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
    
    sys.h2Norm=nrm;
    if inputname(1)
        assignin('caller', inputname(1), sys);
    end
else
    error(['H_' num2str(p) '-norm not implemented.'])
end
end

function [ H_Infty,freq ] = H_Infty( sys )

[mag,w] = freqresp( sys );
[magExtreme]=freqresp(sys,[0,inf]);
mag=cat(3,magExtreme(:,:,1),mag,magExtreme(:,:,2));
w=[0;w;inf];
indices=1:numel(w);
%%Defining Boundaries
maxNorm=sqrt(sum(sum(abs(mag).^2,1),2)); %Frobenius Norm immer smaller than the real norm
[~,indexMaxNorm]=max(maxNorm);
index=indexMaxNorm;
H_Infty=norm(mag(:,:,index),2);
i=0;
%Find the maximum norm computed in freqresp
while(1)
    i=i+1;
    maxNorm(index)=0;
    indices(maxNorm<=H_Infty)=[];
    maxNorm(maxNorm<=H_Infty)=[];
    if not(numel(maxNorm))
        break;
    end
    [~,index]=max(maxNorm);
    computedNorm=norm(mag(:,:,indices(index)),2);
    if (computedNorm>H_Infty)
        H_Infty=computedNorm;
        indexMaxNorm=indices(index);
    end
end
freq=w(indexMaxNorm);

%Optimization of the maximum Hinfty norm
[A,B,C,D,E]=dssdata(sys);
minusA=-A;
w=freq;
[eigenvectors,eigenvalues]=eig(mag(:,:,indexMaxNorm)'*mag(:,:,indexMaxNorm));
[~,Index]=max(diag(eigenvalues));
v=eigenvectors(:,Index);
[Deriv0,Deriv1,Deriv2]=computeDerivatives(minusA,B,C,D,E,w);
lambida=norm(v'*Deriv0)/norm(v');
Vals=[v;lambida;w];
firstDeriv=[-v'*(Deriv0+transp(Deriv0))+2*lambida*v',(v'*v)-1,real(-v'*Deriv1*v)]';
secondDeriv=[-(Deriv0+transp(Deriv0))+2*lambida*eye(size(Deriv0,1)),2*v,(-v'*(Deriv1+transp(Deriv1)))';...
    2*v',0,0;...
    (-v'*(Deriv1+transp(Deriv1))),0,real(-v'*Deriv2*v)];
delta=inf;
i=0;
while(1)
    i=i+1;
    deltaBefore=delta;
    delta=secondDeriv\firstDeriv;
    Vals=Vals-delta;
    w=real(Vals(end));
    lambida=Vals(end-1);
    v=Vals(1:end-2);
    [Deriv0,Deriv1,Deriv2]=computeDerivatives(minusA,B,C,D,E,w);
    firstDeriv=[-v'*(Deriv0+transp(Deriv0))+2*lambida*v',(v'*v)-1,real(-v'*Deriv1*v)]';
    if (abs((norm(deltaBefore)-norm(delta))/norm(delta))<1e-9)|norm(delta(end))<eps(w)|i>=10
        break;
    end
    secondDeriv=[-(Deriv0+transp(Deriv0))+2*lambida*eye(size(Deriv0,1)),2*v,(-v'*(Deriv1+transp(Deriv1)))';...
        2*v',0,0;...
        (-v'*(Deriv1+transp(Deriv1))),0,real(-v'*Deriv2*v)];
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
Deriv1=((respp)'*resp+(resp)'*respp);
%Compute second derivative
b=S\(E*LinearSolve1); b=b(k,:);
LinearSolve2=L\b;
LinearSolve2(l,:)=U\LinearSolve2;
resppp=-2*C*LinearSolve2;
Deriv2=((resppp)'*resp+resp'*resppp+2*(respp)'*respp);
end

