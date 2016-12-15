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
%       If the option 'method' is set to 'adi',then a low rank approximation 
%       of the Cholesky factor [1] is performed. If this option is 
%       specified as 'auto', then ADI is applied to systems with n>500. The  
%       options 'lse', 'rctol' and 'q' only apply to ADI.
%
%       If the option 'rctol' is set, the resulting order can be smaller than
%       specified by the option 'q'. If the option 'forceOrder' is true, the
%       'rctol' is ignored and the result is of order 'q' (not recommended).
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%		*Optional Input Arguments:*
%       -Opts:              a structure containing following options
%           -.method:       select solver for lyapunov equation 
%                           [{'auto'} / 'adi' / 'hammarling' ]
%           -.lse:          solve linear system of equations (only ADI)
%                           [{'gauss'} / 'luChol']
%           -.rctol:        tolerance for difference between ADI iterates
%                           [{'1e-12'} / positive float]
%           -.q:            size of Cholesky factor (only ADI)
%                           [{'0'} / positive integer]
%           -.forceOrder:   return order q
%                           [{'false'} / 'true']
%           -.maxiter:      maximum number of iterations (only ADI)
%                           [{300} / positive integer]
%
% Output Arguments:
%       -R:     Cholesky factor X=R'*R of lyapunov equation A*X+X*A'+B*B'=0
%       -L:     Cholesky factor X=L'*L of lyapunov equation A'*X+X*A+C'*C=0
%
% Examples:
%       Compute the Cholesky factors for both Lyapunov equations
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
%       * *[1] Penzl (2000)*, LYAPACK - A MATLAB Toolbox for Large Lyapunov
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
% Last Change:  11 Dec 2016
% Copyright (c) 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Option parsing
%  Default execution parameters
Def.method  = 'auto';           % 'auto', 'adi', 'hammarling'
Def.lse     = 'gauss';          % only for MESS (see solveLse)
Def.rctol   = 1e-12;            % only for MESS
Def.q       = 0;                % only for MESS
Def.forceOrder  = false;        % only for MESS
Def.maxiter = min([150,sys.n]); % only for MESS

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% check "method" option or make automatic selection

if strcmp(Opts.method,'adi')
    if sys.isDae %DAEs are not supported in sss yet
        warning('Dae-System detected. Using built-in lyapchol instead of ADI.');
        Opts.method='hammarling';
    elseif sys.n<100
        warning('System is small (n<100). Using built-in lyapchol instead of ADI.');
        Opts.method='hammarling';
    end
elseif strcmp(Opts.method,'auto')
    %   Automatic selection of method depending on order and model method
    if sys.n>500 && ~sys.isDae
        Opts.method = 'adi'; %set ADI as default
    else
        Opts.method = 'hammarling';
    end
end

if strcmp(Opts.method, 'adi') % check if MESS is available
    if ~exist('operatormanager.m','file')||~exist('mess_para.m','file')||~exist('mess_lradi.m','file') 
        warning('MESS toolbox not found. Using built-in lyapchol instead of ADI.');
        Opts.method='hammarling';
    end
    pathUsfs=which('operatormanager');
    if strcmp(Opts.method, 'adi') && ~exist(strrep(pathUsfs,'operatormanager.m','solveLse'),'dir')
        warning('MESS user function (usfs) "solveLse" not found. Using usfs "default" instead of "solveLse".');
        Opts.method='hammarling';
        lseType='default';
    else
        lseType='solveLse';
    end
end

%% Solve the lyapunov equation
switch Opts.method
    case 'adi'
    if Opts.q>0 %size of cholesky factor [sys.n x q] -> qmax=q
        Opts.maxiter=Opts.q;
        if Opts.forceOrder
            Opts.rctol=0;
        end
    end   
        
    %% M-MESS ADI
    % eqn struct: system data
    eqn=struct('A_',sys.A,'E_',sys.E,'B',sys.B,'C',sys.C,'prm',speye(size(sys.A)),'type','N','haveE',sys.isDescriptor);
    
    % opts struct: MESS options
    messOpts.adi=struct('shifts',struct('l0',20,'kp',50,'km',25,'b0',ones(sys.n,1),...
        'info',0),'maxiter',Opts.maxiter,'restol',0,'rctol',Opts.rctol,...
        'info',0,'norm','fro');
    
    oper = operatormanager(lseType);
    messOpts.solveLse.lse=Opts.lse;
    messOpts.solveLse.krylov=0;
    
    % get adi shifts
    [messOpts.adi.shifts.p,~,~,~,~,~,~,eqn]=mess_para(eqn,messOpts,oper);
    
    % low rank adi
    [R,Rout,eqn]=mess_lradi(eqn,messOpts,oper);
    
    if Opts.q && size(R,2)<Opts.q
        warning(['Because of small relative changes in the last ADI iterations,',...
            ' the size of R is set to q_R = ',num2str(size(R,2),'%i'),'.']);
    end
    if Rout.rc(end)>Opts.rctol
        warning(['Maximum number of ADI iterations reached (maxiter = ',num2str(Opts.maxiter,'%d'),...
                '). rctol is not satisfied for R: ',num2str(Rout.rc(end),'%d'),' > rctol (',num2str(Opts.rctol,'%d'),').']);
    end
    
    if nargout>1
        if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
            L=R;
        else
            eqn.type='T';
            [L,Lout]=mess_lradi(eqn,messOpts,oper);
        end
        if Opts.q && size(L,2)<Opts.q
            warning(['Because of small relative changes in the last ADI iterations,',...
                ' the size of L is set to q_L = ',num2str(size(L,2),'%i'),'.']);
        end
        if Lout.rc(end)>Opts.rctol
            warning(['Maximum number of ADI iterations reached (maxiter = ',num2str(Opts.maxiter,'%d'),...
                '). rctol is not satisfied for L: ',num2str(Lout.rc(end),'%d'),' > rctol (',num2str(Opts.rctol,'%d'),').']);
        end
    end
    
    % make sure cholesky factorization is like built-in (X=R'*R)
    R=R';
    if nargout>1
        L=L';
    end
    
    case 'hammarling'
        %% built-in lyapchol (hammarling)
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
    otherwise 
        error('sss:lyapchol:invalidMethod','The chosen method for lyapchol is invalid.')
end
end