function [G, omega] = freqresp(varargin)
% FREQRESP - Evaluates complex transfer function of LTI systems
% 
% Description:
%       Evaluates complex transfer function of LTI systems. If the vector
%       of complex frequencies is not passed, then a range of imaginary
%       frequencies is automatically selected.
%
% Syntax:
%       G = freqresp(sys, s)
%       G = freqresp(sys, s, opts)
%       [G, omega] = freqresp(sys)
%
% Inputs:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system or an array with many
%       LTI-Systems with the same number of inputs and outputs.
%       -s: vector of complex frequencies
%       
% Outputs:      
%       -G: vector of complex frequency response values
%       -omega: vector with the frequencies at which the response was computed
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
%       bode, sigma     
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
% Authors:      Jorge Luiz Moreira Silva, Stefan Jaensch, Heiko Panzer, Sylvia Cremer, Rudy Eid
%               Lisa Jeschek, Alessandro Castagnotto
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Dez 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% input parsing
sys= varargin{1};
nOutputs=sys.p;
nInputs=sys.m;
[A,B,C,D,E]=dssdata(sys);
%Reordering Matrices
ACopy=A;
ACopy(find(E))=1;
reOrder=symrcm(ACopy);
sys=sss(A(reOrder,reOrder),B(reOrder,:),C(:,reOrder),D,E(reOrder,reOrder));
%Verifying relation between Inputs and Outputs
M=InputOutputRelation(sys);

omegaIndex = cellfun(@isfloat,varargin);
if ~isempty(omegaIndex) && nnz(omegaIndex)
    %frequency vector was specified
    omega = varargin{omegaIndex};
    varargin(omegaIndex)=[];
else
    %Finding mininum and maximum frequencies
    minW=findminW(sys,M);
    maxW=findmaxW(sys,M);
    %Compute first points of frequency response
    qttyPoints=ceil(log(10^(maxW-minW))/log(2))+1;
    omega=logspace(minW,maxW,qttyPoints)'; %w should be a column according to built-in MATLAB function
    [firstDerivLog,secondDerivLog,magnitude,resp]=ComputeFreqResp(sys,omega*1i,M);
    %Refine the frequency response points
    [G,omega]=FreqRefinement(sys,omega,firstDerivLog,secondDerivLog,magnitude,resp,M);
    return
end

%%  Compute the value of the transfer function at selected freq.
if (sys.Ts==0) % Convert frequency to either laplace or z variable
    s = 1i* omega;
else
    s = exp(1i* omega*sys.Ts);
end

    G=zeros(nOutputs,nInputs*length(varargin),numel(s));
for iSys=1:length(varargin)
    sys= varargin{iSys};
    m=sys.m; p=sys.p;
    if (m~=nInputs|p~=nOutputs)
        error('The number of inputs and outputs of all input-systems must be the same');
    end
    [~,~,~,resp]=ComputeFreqResp(sys,s,M);
    G(1:nOutputs,(1:nInputs)*iSys,:) = resp;
end
end



function [resp,w]=FreqRefinement(sys,w,firstDerivLog,secondDerivLog,magnitude,resp,M)
nInputs=sys.m;
nOutputs=sys.p;
increment=w(2)/w(1);
%Enter a loop for the refinement
while(1)
    %Compute the currentIncrement between all the points
    currentIncrement=reshape(w(2:end)./w(1:end-1),1,1,numel(w)-1);
    currentIncrement=repmat(currentIncrement,nOutputs,nInputs,1);
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
end

function [firstDerivLog,secondDerivLog,magnitude,resp]=ComputeFreqResp(sys,wEval,M)
[A,B,C,D,E]=dssdata(sys);
minusA=-A;
nInputs=sys.m;
nOutputs=sys.p;
resp=zeros(nOutputs,nInputs,numel(wEval));
respp=zeros(nOutputs,nInputs,numel(wEval));
resppp=zeros(nOutputs,nInputs,numel(wEval));
M=repmat(not(M),1,1,numel(wEval));
for i=1:numel(wEval)
    w=wEval(i);
    %Computation of the frequency response for w=wEval(i)
    [L,U,k,l,S]=lu(minusA+E*w,'vector');
    warning('OFF', 'MATLAB:nearlySingularMatrix');
    b=S\B; b=b(k,:);
    linSolve0=L\b;
    linSolve0(l,:)=U\linSolve0;
    resp(:,:,i)=(C*linSolve0)+D;
    %Computation of the first derivative (Respp) of the frequency response for w=wEval(i)
    b=S\(E*linSolve0); b=b(k,:);
    linSolve1=L\b;
    linSolve1(l,:)=U\linSolve1;
    respp(:,:,i)=-w*C*linSolve1;
    %Computation of the second derivative (Resppp) of the frequency response for w=wEval(i)
    b=S\(E*linSolve1); b=b(k,:);
    linSolve2=L\b;
    linSolve2(l,:)=U\linSolve2;
    warning('ON', 'MATLAB:nearlySingularMatrix');
    resppp(:,:,i)=2*C*linSolve2;
end
%Computation of the first two derivatives for a log-log of the magnitude plot
resppp=resppp.*repmat(reshape(wEval,1,1,numel(wEval)),nOutputs,nInputs,1).^2+respp;
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
firstDerivLog(M)=0;
secondDerivLog(M)=0; 
%Warning to guarantee that any Derivative will get to infinity. 
if any(any(any(isinf(firstDerivLog)))) || any(any(any(isinf(secondDerivLog))))
    warning('The magnitude values of your transfer function got too small and the results might not be precise');
end
end

function [M]=InputOutputRelation(sys)
%Function to determine which inputs and outputs are connected through the
%matrices A and E.
%The output is a matrix M p x m.
[A,B,C,~,E]=dssdata(sys);
C(find(C))=1;
A(find(E))=1;
vec=zeros(sys.m,size(A,1));
 for k=1:sys.m
    [i,~]=find(B(:,k));
    vec(k,i)=1;
while(1)
    [~,j]=find(A(i,:));
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
    %Matrix=H+identity*minW*1i;
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
%Gershgorin circle theorem
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
