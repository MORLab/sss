function [varargout] = freqresp(varargin)
% FREQRESP - Frequency response of sparse state-space systems.
% 
% Syntax:
%       [G, omega] = freqresp(sys)
%       G = freqresp(sys, omega)
%       G = freqresp(sys, Opts)
%       G = freqresp(sys, omega, Opts)
%
% Description:
%       Evaluates complex transfer function of LTI systems. 
%
%       If the vector of complex frequencies is not passed, then a range of 
%       imaginary frequencies is automatically selected. 
%       This automatic selection is done through the computation of the
%       first and second derivatives of the magnitude of the frequency
%       response. 
%
% Inputs:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -omega: vector of frequencies or cell with {wmin,wmax}
%       -Opts:  structure with execution parameters
%           -.maxPoints: Maximum number of refinement points
%                       [{1500} / positive integer]
%           -.lse:  solve linear system of equations
%                       [{'sparse'} / 'full' /'gauss' /'hess' / 'iterative']
%       
% Outputs:      
%       -G: |sys.p| x |sys.m| x |N| array of complex frequency response values,
%       where |N| is the number of sampling frequencies
%       -w: vector with the frequencies at which the response was computed
%
% Examples:
%       The following code computes the frequency response of the benchmark
%       'building' and returns also the vector of frequencies at which the
%       response was computed |omega|:
%
%> load building.mat
%> sys=sss(A,B,C);
%> [G,omega]=freqresp(sys);
%
% See Also:
%       bode, sigma, bodemag, bodeplot, frd
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
% Authors:      Jorge Luiz Moreira Silva, Stefan Jaensch, Heiko Panzer, Sylvia Cremer, Rudy Eid
%               Lisa Jeschek, Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  14 Sep 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Parse inputs and options
Def.maxPoints = 1500; % maximum number of refinement points
Def.lse = 'sparse'; % solveLse

% create the options structure
if ~isempty(varargin) && isstruct(varargin{end})
    Opts        = varargin{end};
    varargin    = varargin(1:end-1);
end
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% Frequency vector
omega           = [];
omegaIndex      = cellfun(@isfloat,varargin);
omegaCellIndex  = cellfun(@iscell, varargin);
if ~isempty(omegaIndex) && nnz(omegaIndex)
    omega       = varargin{omegaIndex};
    varargin(omegaIndex)=[];
elseif ~isempty(omegaCellIndex) && nnz(omegaCellIndex)
    minW        =log10(varargin{omegaCellIndex}{1});
    maxW        =log10(varargin{omegaCellIndex}{2});
    varargin(omegaCellIndex)=[];
end

sys= varargin{1};

if sys.isDae == true && isempty(omega)
    error('freqresp and bode only support DAEs, if a single frequency or a frequency vector omega is parsed to the functions');
end

nOutputs        = sys.p;
nInputs         = sys.m;
[A,B,C,D,E]     = dssdata(sys);
%Reordering Matrices
reOrderMatrix   = abs(A)+abs(E); %abs is used to guarantee that any terms will be cancelled in this combination of A and E
reOrder         = symrcm(reOrderMatrix);
% reOrderMatrix   = reOrderMatrix(reOrder,reOrder);
sys.A           = A(reOrder,reOrder);
sys.B           = B(reOrder,:);
sys.C           = C(:,reOrder);
sys.E           = E(reOrder,reOrder);

%% Verifying relation between Inputs and Outputs

%{ 
One of the steps is the computation of the connection between inputs and 
outputs of the system. The value of M(i,j) is one when the input j is 
connected to the output i. When it is zero, then a change in the input j 
doesn't influence the output i.
%}
M=InputOutputRelation(sys);

if ~any(any(M))
    %Disconnected, but static gain
    if not(exist('omega','var')) || isempty(omega)
        if isempty(omegaCellIndex) || ~nnz(omegaCellIndex)
            maxW=1;
            minW=0;
        end
        qttyPoints=ceil(log(10^(maxW-minW))/log(2))+1;
        omega=logspace(minW,maxW,qttyPoints)';
    else
        %make sure it's a column vector
        if size(omega,2)>size(omega,1)
            omega = omega.';
        end
    end
    G=zeros(nOutputs,nInputs,size(omega,1));
    for i=1:nOutputs
        for j=1:nInputs
            G(i,j,:)=ones(size(omega,1),1)*sys.D(i,j);
        end
    end
    
elseif not(exist('omega','var')) || isempty(omega)
    %Finding mininum and maximum frequencies
    if isempty(omegaCellIndex) || ~nnz(omegaCellIndex)
        minW    =findminW(sys,M);
        maxW    =findmaxW(sys,M);
    end
    %Compute first points of frequency response
    if minW == -Inf, minW = findminW(sys,M); end
    qttyPoints=ceil((maxW-minW)/log(2))+1;
    omega = real(logspace(minW,maxW,qttyPoints))'; %omega should be a column according to built-in MATLAB function
    [firstDerivLog,secondDerivLog,magnitude,resp]=ComputeFreqResp(sys,omega*1i,M);
    %Refine the frequency response points
    [G,omega]   = FreqRefinement(sys,omega,firstDerivLog,secondDerivLog,magnitude,resp,M,Opts);
else
    %%  Compute the value of the transfer function at selected freq.
    %make sure it's a column vector
    if size(omega,2)>size(omega,1), omega = omega.'; end
    
    if (sys.Ts==0) % Convert frequency to either laplace or z variable
        %{
        % if omega is a real vector, then make it imaginary for the transfer
        % function evaluations. If it is purely imaginary, complex or mixed,
        % than just use the values passed.
        % NOTE: this is in line with the built-in case. However, if one desired
        % to evaluate the transfer function ONLY on the real axis, this implementation does not allow it. 
        % In this case, we recommend evaluating the transfer function manually
        % or add an imaginary frequency which can be disregarded afterwards
        %}
        if isreal(omega)
            s = 1i* omega;
        else
            s = omega;
        end
    else
        if isreal(omega)
            s = exp(1i* omega*sys.Ts);
        else
            s = omega;
        end
    end

    G=zeros(nOutputs,nInputs,numel(s));
    [~,~,~,resp]=ComputeFreqResp(sys,s,M);
    G(1:nOutputs,1:nInputs,:) = resp;
end

varargout{1}=G;
varargout{2}=omega;


end



function [resp,w]=FreqRefinement(sys,w,firstDerivLog,secondDerivLog,magnitude,resp,M,Opts)
nInputs=sys.m;
nOutputs=sys.p;
increment=w(2)/w(1);
%Enter a loop for the refinement
while(length(w)<=Opts.maxPoints)
    %Compute the currentIncrement between all the points
    currentIncrement=reshape(w(2:end)./w(1:end-1),1,1,numel(w)-1);
    currentIncrement=repmat(currentIncrement,nOutputs,nInputs);
    %Compute expected value considering just first derivative (make the
    %plot smooth).
    Expected=exp(log(magnitude(:,:,1:end-1))+log(currentIncrement).*firstDerivLog(:,:,1:end-1));
    delta1=abs((Expected-magnitude(:,:,2:end))./(magnitude(:,:,2:end)));
    delta1(isnan(delta1))=0;
    %Compute expected value considering first and second derivatives
    Expected=exp(log(magnitude(:,:,1:end-1))+log(currentIncrement).*firstDerivLog(:,:,1:end-1)+log(currentIncrement).^2.*secondDerivLog(:,:,1:end-1)/2);
    delta2=abs((Expected-magnitude(:,:,2:end))./(magnitude(:,:,2:end)));
    delta2(isnan(delta2))=0;
    %Conditions to refine. If the values get smaller, the final bode will
    %be more refined.
    wEvaluated=w(any(any(or(delta1>0.05,delta2>0.01),1),2));
    if not(numel(wEvaluated))
        break;
    end
    %Compute new incrmeent
    increment=sqrt(increment);
    if (increment==1) %it is not possible to increment anymore.
        break;
    end
    wEvaluated=wEvaluated*increment;
    %Compute the new frequency responses
    [firstDerivLogComp,secondDerivLogComp,magnitudeComp,respComp]=ComputeFreqResp(sys,wEvaluated*1i,M);
    %Reorder everything
    [w,Index]=sort([w;wEvaluated]);
    firstDerivLog=cat(3,firstDerivLog,firstDerivLogComp);
    firstDerivLog=firstDerivLog(:,:,Index);
    secondDerivLog=cat(3,secondDerivLog,secondDerivLogComp);
    secondDerivLog=secondDerivLog(:,:,Index);
    magnitude=cat(3,magnitude,magnitudeComp);
    magnitude=magnitude(:,:,Index);
    resp=cat(3,resp,respComp);
    resp=resp(:,:,Index);
end
if length(w)>=Opts.maxPoints
    warning(['Maximum number of refinement points reached. Increase Opts.maxPoints '...
        'for a better resolution.']);
end
end

function [firstDerivLog,secondDerivLog,magnitude,resp]=ComputeFreqResp(sys,wEval,M)
[A,B,C,D,E]=dssdata(sys);
minusA=-A;
nInputs=sys.m;
nOutputs=sys.p;
resp=zeros(nOutputs,nInputs,numel(wEval));
respp=zeros(nOutputs,nInputs,numel(wEval));
resppp=zeros(nOutputs,nInputs,numel(wEval));
if length(M)==1 && nnz(M)
    M=zeros(1,1,numel(wEval));
else
    M=repmat(not(M),1,1,numel(wEval));
end
for i=1:numel(wEval)
    w=wEval(i);
    if isinf(w)
        resp(:,:,i)=D;
        respp(:,:,i)=nan(size(D));
        resppp(:,:,i)=nan(size(D));
    else
        %Computation of the frequency response for w=wEval(i)
        Opts.krylov='standardKrylov';

        % tangential directions
        Rt=zeros(size(B,2),size(B,2)*3);
        for nB2=1:size(B,2)
            Rt(nB2,(nB2-1)*3+1:(nB2-1)*3+3)=ones(1,3);
        end
        
        % solve lse
        linSolve=solveLse(minusA,B,E,-ones(1,3*size(B,2))*w,Rt,Opts);

        % sort solutions
        linSolve1=zeros(size(B));
        linSolve2=zeros(size(B));
        linSolve3=zeros(size(B));
        for nB2=1:size(B,2)
            linSolve1(:,nB2)=linSolve(:,(nB2-1)*3+1);
            linSolve2(:,nB2)=linSolve(:,(nB2-1)*3+2);
            linSolve3(:,nB2)=linSolve(:,(nB2-1)*3+3);
        end
        
        resp(:,:,i)=(C*linSolve1)+D;
        
        %Computation of the first derivative (Respp) of the frequency response for w=wEval(i)
        respp(:,:,i)=-w*C*linSolve2;

        %Computation of the second derivative (Resppp) of the frequency response for w=wEval(i)
        resppp(:,:,i)=2*C*linSolve3;
    end
end
%Computation of the first two derivatives for a log-log of the magnitude plot
resppp=resppp.*repmat(reshape(wEval,1,1,numel(wEval)),nOutputs,nInputs).^2+respp;
limit=sqrt(eps(0))*10^5; %Limit for change the order of computation (see below)
magnitude=abs(resp);
cond=(any(any(magnitude<limit,1),2));
magSquared=zeros(nOutputs,nInputs,numel(wEval));
firstDerivLog=zeros(nOutputs,nInputs,numel(wEval));
secondDerivLog=zeros(nOutputs,nInputs,numel(wEval));
%Computation of first and second derivatives of conj(Resp).*Resp
Deriv1=(conj(respp).*resp+conj(resp).*respp);
Deriv2=2*(real(conj(resppp).*resp)+(conj(respp).*respp));

i=0; %Computation of log derivatives considering only magnitudes>=limit
w=(cond==i);
magSquared(:,:,w)=(conj(resp(:,:,w)).*resp(:,:,w));
firstDerivLog(:,:,w)=(0.5*Deriv1(:,:,w)./magSquared(:,:,w));
secondDerivLog(:,:,w)=0.5*(Deriv2(:,:,w).*magSquared(:,:,w)-Deriv1(:,:,w).^2)./magSquared(:,:,w)./magSquared(:,:,w); 

i=1; %Computation of log derivatives considering only magnitude<limit
w=(cond==i);
firstDerivLog(:,:,w)=(0.5*Deriv1(:,:,w)./magnitude(:,:,w)./magnitude(:,:,w));
secondDerivLog(:,:,w)=0.5*(Deriv2(:,:,w).*magnitude(:,:,w).*magnitude(:,:,w)-Deriv1(:,:,w).^2)./magnitude(:,:,w)./magnitude(:,:,w)./magnitude(:,:,w)./magnitude(:,:,w); 

%The derivatives of transfer functions of inputs and outputs not connected must be zero
if nnz(M)
    firstDerivLog(M)=0;
    secondDerivLog(M)=0; 
end
%Warning to guarantee that any Derivative will get to infinity. 
if any(any(any(isinf(firstDerivLog)))) || any(any(any(isinf(secondDerivLog))))
    warning('The magnitude values of your transfer function got too small and the results might not be precise');
end
end

function [M]=InputOutputRelation(sys)
%Function to determine which inputs and outputs are connected through the
%matrices A and E.
%The output is a matrix M p x m.
[~,B,C,~,~]=dssdata(sys);
C(find(C))=1;
% vec=zeros(sys.m,size(S,1));
vec=zeros(sys.m,sys.n);
for k=1:sys.m
    [i,~]=find(B(:,k));
    vec(k,i)=1;
    while(1)
%         [j,~]=find(S(:,i));
        [~,K]=find(sys.E(i,:));
        [j,~]=find(sys.A(:,K));
        j=unique(j);
        newVec=(not(vec(k,j)));
        vec(k,j(newVec))=1;
        i=j(newVec);
        if numel(newVec)==0
            break;
        end
    end
end
M=C*vec';
M(find(M))=1;
end

function [minWFinal] = findminW(sys,M)
[A,B,C,~,E]=dssdata(sys);
minusA=-A;
nInputs=size(B,2);
nOutputs=size(C,1);

normE=norm(E,1);
normInvA=condest(minusA)/norm(minusA,1);
minW=1/(normInvA*normE);
if minW>100*eps
    minW=minW/100;
end
% finding minW
if minW<=100*eps
    minW=100*eps;
end

i=0;
minW=minW/10;
CondSuf=not(M);
Cond=zeros(1,2);
%minW
secondDerivLog=0;
while(1)
    i=i+1;
    if i>1
        minW=minW*10;
    end
    Matrix=minusA+E*minW*1i;
    if any(any(CondSuf==0))
        Cond=repmat(condest(Matrix),nOutputs,nInputs);
    end
    secondDerivLogBefore=secondDerivLog;
    [~,secondDerivLog,~,~]=ComputeFreqResp(sys,minW*1i,M);
    Var=inf(nOutputs,nInputs);
    Var(M&(CondSuf|(Cond>10^-3)))=1e-3;
    if i>1
        if Cond>100*eps
            CondSuf((abs(secondDerivLog)>abs(secondDerivLogBefore))&not(CondSuf))=1;
        end
        if any(any(abs(secondDerivLog)>Var))
            break;
        end
    end
    
end
minWFinal=floor(log10(minW(end)*10));

end

function [maxWFinal] = findmaxW(sys,M)
[A,B,C,~,E]=dssdata(sys);
minusA=-A;
nInputs=size(B,2);
nOutputs=size(C,1);
normInvE=condest(E)/norm(E,1);
%Gershgorin circle theorem: one of its implications is that the first norm
%of a matrix will be greater or equal to its biggest absolute eigenvalue.
normA=norm(minusA,1);
maxW=normA*normInvE;

maxW=maxW*10^4;

i=0;
CondSuf=not(M);%zeros(nInputs,nOutputs);
Cond=zeros(1,2);
secondDerivLog=0;
while(1)
    i=i+1;
    if i>1
        maxW=maxW/10;
    end
    Matrix=minusA+E*maxW*1i;
    if any(any(CondSuf==0))
        Cond=repmat(condest(Matrix),nOutputs,nInputs);
    end
    secondDerivLogBefore=secondDerivLog;
    [~,secondDerivLog,~,~]=ComputeFreqResp(sys,maxW*1i,M);
    Var=inf(nOutputs,nInputs);
    Var(M&(CondSuf|(Cond>10^-3)))=1e-3;
    if i>1
        if Cond>100*eps
            CondSuf((abs(secondDerivLog)>abs(secondDerivLogBefore))&not(CondSuf))=1;
        end
        if any(any(abs(secondDerivLog)>Var))
            break;
        end
    end
    
end
maxWFinal=floor(log10(maxW(end)*10));
end
