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
%       -p: choice of H_2-norm or H_inf-norm 
%           [{'2'} / 'inf']
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
% Authors:      Heiko Panzer, Sylvia Cremer, Rudy Eid
%               Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Feb 2016
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
        [sys.hInfNorm, sys.hInfPeakfreq] = norm(ss(sys),inf);
        if inputname(1)
            assignin('caller', inputname(1), sys);
        end
    end
    nrm=sys.hInfNorm; 
    if nargout>1
        varargout{1}=sys.hInfPeakfreq;
    end
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
                    if ~sys.isDae
                        % Lyapack opts for Adi
                        lyaOpts.l0=20;
                        lyaOpts.kp=50;
                        lyaOpts.km=25;
                        lyaOpts.method='heur';
                        lyaOpts.zk='Z';
                        lyaOpts.rc='C';
                        lyaOpts.adi=struct('type','B','max_it', 100,'min_res',0,'with_rs','N',...
                            'min_in',1e-12,'min_end',0,'info',0,'cc_upd',0,'cc_tol',0);
                    end
                    if sys.isDescriptor
                        try
                            if sys.n<100
                                error('System is too small for ADI. ');
                            end
                            if sys.isDae
                                error('ADI does not work with DAE systems. ');
                            end
                            if isstable(sys)~=1
                                warning('System appears to be unstable. The norm will be set to Inf.');
                                nrm=Inf;
                                return;
                            end
                            if sys.isSym
                                 lyaOpts.usfs=struct('s','msns_s','m','msns_m');
                                [M0,MU0,N0,B0,C0]=msns_pre(sys.E,sys.A,sys.B,sys.C);
                                msns_m_i(M0,MU0,N0); 
                                msns_l_i;
                                p=lp_para(msns,[],[],lyaOpts,ones(length(B0),1));
                                lyaOpts.p=p.p;
                                msns_s_i(lyaOpts.p);
                            else
                                lyaOpts.usfs=struct('s','munu_s','m','munu_m');
                                [M0,ML0,MU0,N0,B0,C0]=munu_pre(sys.E,sys.A,sys.B,sys.C);
                                munu_m_i(M0,ML0,MU0,N0); 
                                munu_l_i;
                                p=lp_para(munu,[],[],lyaOpts,ones(length(B0),1));
                                lyaOpts.p=p.p;
                                munu_s_i(lyaOpts.p);
                            end

                            R=lp_lradi([],[],B0,lyaOpts); %ADI solution of lyapunov equation
                            nrm=norm(R'*C0','fro');
                            munu_m_d;
                            munu_l_d;
                            munu_s_d(p.p);
                        catch ex
                            warning([ex.message,'Trying without ADI...']);
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
                        end
                    else
                        try
                            if sys.n<100
                                error('System is too small for ADI. ');
                            end
                            if isstable(sys)~=1
                                warning('System appears to be unstable. The norm will be set to Inf.');
                                nrm=Inf;
                                return;
                            end
                            if sys.isSym
                                lyaOpts.usfs=struct('s','as_s','m','as_m');
                                [A0,B0,C0]=as_pre(sys.A,sys.B,sys.C);
                                as_m_i(A0);
                                as_l_i;
                                p=lp_para(as,[],[],lyaOpts, ones(length(B0),1));
                                lyaOpts.p=p.p;
                                as_s_i(lyaOpts.p);
                            else
                                lyaOpts.usfs=struct('s','au_s','m','au_m');
                                [A0,B0,C0]=au_pre(sys.A,sys.B,sys.C);
                                au_m_i(A0);
                                au_l_i;
                                p=lp_para(au,[],[],lyaOpts, ones(length(B0),1));
                                lyaOpts.p=p.p;
                                au_s_i(lyaOpts.p);
                            end

                            R=lp_lradi([],[],B0,lyaOpts); %ADI solution of lyapunov equation
                            nrm=norm(R'*C0','fro');
                            au_m_d;
                            au_l_d;
                            au_s_d(p.p);
                        catch ex
                            warning([ex.message,'Trying without ADI...']);
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
    
    sys.h2Norm=nrm;
    if inputname(1)
        assignin('caller', inputname(1), sys);
    end
else
    error(['H_' num2str(p) '-norm not implemented.'])
end
