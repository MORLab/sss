function [R,L] = lyapchol(sys, Opts)
% LYAPCHOL - Solve Lyapunov equations
%
% Syntax:
%       R				= LYAPCHOL(sys)
%       [R,L]			= LYAPCHOL(sys)
%       [R,L]  		    = LYAPCHOL(sys,Opts)
%
% Description:
%       This function returns the Cholesky factorization X=R'*R of the 
%       solution of the Lyapunov equation A*X+X*A'+B*B'=0.
%
%       If the option 'type' is set to 'adi',then a low rank approximation 
%       of the Cholesky factor [1] is performed. If this option is not 
%       specified, then ADI is applied to systems with n>500. The options 
%       'lse', 'rctol' and 'q' only apply to ADI.
%
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%		*Optional Input Arguments:*
%       -Opts:              a structure containing following options
%           -.type:         select amongst different tbr algorithms
%                           [{''} / 'adi' / 'builtIn' ]
%           -.lse:          solve linear system of equations (only ADI)
%                           [{'gauss'} / 'luChol']
%           -.rctol:        tolerance for difference between ADI iterates
%                           [{'1e-12'} / positive float]
%           -.q:            size of Cholesky factor (only ADI)
%                           [{'0'} / positive integer]
%
% Output Arguments:
%       -R:     Cholesky factor X=R'*R of lyapunov equation A*X+X*A'+B*B'=0
%       -L:     Cholesky factor X=L'*L of lyapunov equation A'*X+X*A+C'*C=0
%
% Examples:
%       Compute the Cholesky factors for both lyapunov equations
%
%> sys = loadSss('building');
%> [R,L] = lyapchol(sys);
%
%       To compute a single Cholesky factor, use
%
%> R = lyapchol(sys);
%> L = lyapchol(sys');
%
% See Also:
%       solveLse, tbr, norm, isrk
%
% References:
%       * *[1] Penzl (2000)*, |lyapack| - A MATLAB Toolbox for Large Lyapunov
%       and Riccati Equations, Model Reduction Problems, and Linear-Quadratic 
%       Optimal Control Problems.
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
% Authors:      Alessandro Castagnotto, 
%               Lisa Jeschek
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Aug 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%   Default execution parameters
Def.type    = ''; % 'adi', 'builtIn'
Def.lse     = 'gauss'; % only for mess ('gauss', 'luChol')
Def.rctol   = 1e-12; % only for mess
Def.q       = 0; % only for mess

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

if strcmp(Opts.type,'adi') && sys.isDae
    warning('MESS does not support dae-Systems. The built-in lyapchol is used instead.');
    Opts.type='';
end

if strcmp(Opts.type,'adi') && sys.n<100
    warning('System is too small for the use of ADI, the built-in lyapchol is used instead.');
    Opts.type='';
end

if (strcmp(Opts.type,'') && sys.n>500 && ~sys.isDae) || strcmp(Opts.type,'adi')
    %% mess
    % eqn struct: system data
    eqn=struct('A_',sys.A,'E_',sys.E,'B',sys.B,'C',sys.C,'prm',speye(size(sys.A)),'type','N','haveE',sys.isDescriptor);
    
    % opts struct: mess options
    messOpts.adi=struct('shifts',struct('l0',20,'kp',50,'km',25,'b0',ones(sys.n,1),...
        'info',0),'maxiter',300,'restol',0,'rctol',1e-12,...
        'info',0,'norm','fro');
    
    oper = operatormanager('solveLse');
    messOpts.solveLse.lse=Opts.lse;
    messOpts.solveLse.krylov=0;
    
    if Opts.q>0 %size of cholesky factor [sys.n x q] -> qmax=q
        messOpts.adi.maxiter=Opts.q;
        messOpts.adi.restol=0;
        messOpts.adi.rctol=1e-30;
    end
    
    % get adi shifts
    [messOpts.adi.shifts.p, eqn]=mess_para(eqn,messOpts,oper);
    
    % low rank adi
    [R,Rout,eqn]=mess_lradi(eqn,messOpts,oper);
    
    if nargout>1
        if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
            L=R;
        else
            eqn.type='T';
            [L,Lout]=mess_lradi(eqn,messOpts,oper);
        end
    end
    
    if Opts.q>0 % warn user if rctol is satisfied before q_user
        qminR=Opts.q;
        qminL=Opts.q;
        nStop=0;

        % rctol is satisfied if rc<tol for 10 times consecutively
        for i=1:length(Rout.rc)
            if Rout.rc(i)<Opts.rctol
                nStop=nStop+1;
            else
                nStop=0;
            end
            if nStop==10
                qminR=i;
                break
            end
        end

        if nargout>1
            if ~(sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0) && qminR<Opts.q
                qminL=Opts.q;
                nStop=0;
                for i=1:length(Lout.rc)
                    if Lout.rc(i)<Opts.rctol
                        nStop=nStop+1;
                    else
                        nStop=0;
                    end
                    if nStop==10
                        qminL=i;
                        break
                    end
                end
            elseif sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                qminL=qminR;
            end
        else
            qminL=qminR;
        end
        q_min_in=max(qminR,qminL);

        if q_min_in>0 && q_min_in < Opts.q
            warning(['After q=', num2str(q_min_in,'%d'),...
            ' the contribution of the ADI iterates was very small. Consider reducing the desired order accordingly.']);
        end
    end
    
    % make sure cholesky factorization is like built-in (X=R'*R)
    R=R';
    if nargout>1
        L=L';
    end
else
    %% built-in
    if sys.isDescriptor
        R = lyapchol(sys.A,sys.B,sys.E);
    else
        R = lyapchol(sys.A,sys.B);
    end

    if nargout>1
        if sys.isDescriptor
            L = lyapchol(sys.A', sys.C',sys.E');
        else
            L = lyapchol(sys.A',sys.C');
        end
    end    
end
end