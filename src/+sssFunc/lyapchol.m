function [S,varargout] = lyapchol(varargin)
% LYAPCHOL - Solve Lyapunov equations
%
% Syntax:
%       [S,data]        = LYAPCHOL(sys)
%       [S,R,data]		= LYAPCHOL(sys)
%       S               = LYAPCHOL(A,B)
%       S               = LYAPCHOL(A,B,E)
%       ...             = LYAPCHOL(...,Opts)
%
% Description:
%       This function returns the Cholesky factorization X=S*S' of the 
%       solution of the Lyapunov equation A*X+X*A'+B*B'=0 or the generalized 
%       Lyapunov equation A*X*E'+E*X*A'+B*B'=0.
%
%       If the number of output arguments is 2, then the low rank factor 
%       Y = R*R' of the dual (generalized) lyapunov equation 
%       A'*Y*E+E'*Y*A+C'*C=0 is solved as well.
%
%       To call this version of lyapchol with matrices A,B,E, make sure to
%       use |sssFunc.lyapchol(...)|.
%
%       If the option 'type' is set to 'adi',then a low rank approximation 
%       of the Cholesky (like) factor is performed [1]. If this option is not 
%       specified, then ADI is applied to systems with n>500. The options 
%       'lse', 'rctol' and 'q' only apply to ADI.
%
%       //Note: the definition of the Cholesky factors X = S*S' is
%       different from built-in lyapchol, where X = S'*S. However, our
%       definition is consistent both with standard literature (cp [3]) and
%       the low-rank approximation if R has fewer columns than rows.
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%       -A,B,E: (alternatively) matrices of the Lyapunov equation
%		*Optional Input Arguments:*
%       -Opts:                  a structure containing following options
%           -.method:           select solver for lyapunov equation 
%                               [{'auto'} / 'adi' / 'hammarling' / 'crksm' ]
%           -.lse:              solve linear system of equations (only ADI and CRKSM)
%                               [{'gauss'} / 'luChol']
%           -.rctol:            tolerance for difference between ADI or crksm iterates
%                               [{'1e-12'} / positive float]
%           -.restol:           tolerance for the residual of the Lyapunov eqution
%                               [{'1e-8'} / positive float]
%           -.maxiter:          maximum number of iterations (only ADI and CRKSM)
%                               [{150} / positive integer]
%           -.adiShiftsMethod:  refer to opts.adi.shifts.method in
%                               mess_para or mess_lradi
%           -.q:                size of Cholesky factor (only ADI)
%                               [{'0'} / positive integer]
%           -.forceOrder:       return order q
%                               [{'false'} / 'true']
%           -.subspace:         build block or tangential Krylov subspace (only CRKSM)
%                               [{'block'} / 'tangential']
%           -.projection:       choose projection method to get the reduced system (only CRKSM)
%                               [{'onesided'} / 'twosided']
%           -.initShifts:       choose initial shift strategy in INITIALIZESHIFTS (only CRKSM)
%                               [{'ADI'} / 'eigs' / 'ROM' / 'const']
%           -.nShifts:          set number of initial shifts (only CRKSM)
%                               positive, even integer
%           -.shifts:           choose, wether the shifts should be updatet
%                               within the iterations (dynamical) or uses use
%                               the initial set again during  the whole process
%                               (fixedCyclic)
%                               [{'dynamical'} / 'fixedCyclic']
%
% Output Arguments:
%       -S:     Cholesky factor X=S*S' of Lyapunov equation A*X*E'+E*X*A'+B*B'=0
%       -R:     Cholesky factor Y=R*R' of Lyapunov equation A'*Y*E+E'*Y*A+C'*C=0
%
% Examples:
%       Compute the Cholesky factors for both Lyapunov equations
%
%> sys = loadSss('building');
%> [S,R] = lyapchol(sys);
%
%       To compute a single Cholesky factor, use
%
%> S = lyapchol(sys);
%> R = lyapchol(sys');
%
%       To call this version of lyapchol with matrices use
%
%> [A,B,C,D,E] = dssdata(sys);
%> S = sssFunc.lyapchol(A,B,E);
%
% See Also:
%       solveLse, tbr, norm, numerics/lyapchol
%
% References:
%       * *[1] Penzl (2000)*, LYAPACK - A MATLAB Toolbox for Large Lyapunov
%       and Riccati Equations, Model Reduction Problems, and Linear-Quadratic 
%       Optimal Control Problems.
%       * *[2] Saak, Köhler, Benner (2016)*, M-M.E.S.S. - The Matrix 
%       Equation Sparse Solver Library.
%       * *[3] Golub, Van Loan (1996)*, Matrix computations
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
% Authors:      Alessandro Castagnotto, Lisa Jeschek, Paul Heidenreich,
%               Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  21 Nov 2017
% Copyright (c) 2016,2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Input parsing
% Opts
if ~isempty(varargin) && isstruct(varargin{end})
    Opts = varargin{end};
    varargin = varargin(1:end-1);
else
    Opts = struct();
end

% System/Matrices
if isa(varargin{1},'sss') || isa(varargin{1},'ss')
    sys = varargin{1};
else
    %create a mock system to perform computations
    A = varargin{1};
    B = varargin{2};
    if length(varargin)>2
        E = varargin{3};
    else
        E = speye(size(A));
    end
    sys = sss(A,B,sparse(1,size(A,2)),[],E);
end 
%% Option parsing
% General default execution parameters
Def.method          = 'auto';           % 'auto', 'adi', 'hammarling', 'crksm'
Def.lse             = 'gauss';          % only for MESS (see solveLse) and CRKSM
Def.restol          = 1e-8;             % only for MESS and CRKSM
Def.rctol           = 0;                % only for MESS and CRKSM
Def.maxiter         = min([150,sys.n]); % only for MESS and CRKSM
Def.infoLyap        = 0;                % give output data-struct: [{0}, 1]

% ADI default execution parameters
Def.adiShiftsMethod = 'projection';     % only for MESS  
Def.q               = 0;                % only for MESS
Def.forceOrder      = false;            % only for MESS

% CRKSM default execution parameters
Def.subspace        = 'block';          % only for CRKSM; build block or tangential Krylov subspace, [{'block'} / 'tangential']
Def.initShifts      = 'ADI';            % only for CRKSM; choose shift strategy in INITIALIZESHIFTS, [{'ADI'} / 'eigs' / 'ROM' / 'const']
Def.nShifts         = 10;               % only for CRKSM; set number of initial shifts
Def.shifts          = 'dynamical';      % only for CRKSM; [{'dynamical'} / 'fixedCyclic']

% create the options structure
if ~exist('Opts','var') || isempty(fieldnames(Opts))
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% check "method" option or make automatic selection

if strcmp(Opts.method,'adi') ||  strcmp(Opts.method,'crksm')
    if sys.isDae %DAEs are not supported in sss yet
        warning('Dae-System detected. Using built-in lyapchol instead of ADI/CRKSM.');
        Opts.method='hammarling';
    elseif sys.n<100
        warning('System is small (n<100). Using built-in lyapchol instead of ADI/CRKSM.');
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

 % check if CRKSM and internal functions are available
if strcmp(Opts.method, 'crksm') && (~exist('crksm.m','file') || ~exist('initializeShifts.m','file') || ~exist('getShifts.m','file')...
   || ~exist('solveLse.m','file')) 
    warning('One of the following files not found.');
    warning('crksm.m, initializeShifts.m, getShifts.m, solveLse.m');
    warning('Using ADI method instead of CRKSM.');
    Opts.method='adi'; 
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
        if Opts.forceOrder, Opts.rctol=0; end  
        
        %% M-MESS ADI
        % eqn struct: system data
        eqn=struct('A_',sys.A,'E_',sys.E,'B',sys.B,'C',sys.C,'prm',speye(size(sys.A)),'type','N','haveE',sys.isDescriptor);

        % opts struct: MESS options
        messOpts.adi=struct('shifts',struct('l0',20,'kp',50,'km',25,'b0',ones(sys.n,1),'method',Opts.adiShiftsMethod,...
            'info',0),'maxiter',Opts.maxiter,'restol',Opts.restol,'rctol',Opts.rctol,...
            'info',0,'norm','fro');

        oper = operatormanager(lseType);
        messOpts.solveLse.lse=Opts.lse;
        messOpts.solveLse.krylov=0;

        % get adi shifts
        [messOpts.adi.shifts.p]=mess_para(eqn,messOpts,oper);

        % low rank adi
        %messOpts.adi.shifts.method = 'projection';
        [S,Sout]=mess_lradi(eqn,messOpts,oper);
        data.Info_S = Sout;

        if Opts.q && size(S,2)<Opts.q
            warning(['Because of small relative changes in the last ADI iterations,',...
                ' the size of S is set to q_S = ',num2str(size(S,2),'%i'),'.']);
        end
        
        if Sout.niter >= Opts.maxiter
            warning(['Maximum number of ADI iterations reached (maxiter = ',num2str(Opts.maxiter,'%d'), ').']);
        elseif isfield(Sout,'res') && Sout.res(end)>Opts.restol
            warning(['restol is not satisfied for S: ',num2str(Sout.res(end),'%d'),' > rctol (',num2str(Opts.restol,'%d'),').']);
        elseif isfield(Sout,'rc') && Sout.rc(end)>Opts.rctol
             warning(['rctol is not satisfied for S: ',num2str(Sout.rc(end),'%d'),' > rctol (',num2str(Opts.rctol,'%d'),').']);
        end
        
        if (nargout>1 && Opts.infoLyap == 1) || (nargout == 2 && Opts.infoLyap == 0) || nargout == 3
            if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                R=S;
                data.Info_R = Sout;
            else
                eqn.type='T';
                [R,Rout]=mess_lradi(eqn,messOpts,oper);

                data.Info_R = Rout;
                if Rout.niter >= Opts.maxiter
                    warning(['Maximum number of ADI iterations reached (maxiter = ',num2str(Opts.maxiter,'%d'), ').']);
                elseif isfield(Rout,'res') && Rout.res(end)>Opts.restol
                    warning(['restol is not satisfied for R: ',num2str(Rout.res(end),'%d'),' > rctol (',num2str(Opts.restol,'%d'),').']);
                elseif isfield(Rout,'rc') && Rout.rc(end)>Opts.rctol
                    warning(['rctol is not satisfied for R: ',num2str(Rout.rc(end),'%d'),' > rctol (',num2str(Opts.rctol,'%d'),').']);
                end
            end
            if Opts.q && size(R,2)<Opts.q
                warning(['Because of small relative changes in the last ADI iterations,',...
                ' the size of R is set to q_R = ',num2str(size(R,2),'%i'),'.']);
            end
        end
    
    case 'hammarling'
        %% built-in lyapchol (hammarling)
        if sys.isDescriptor
            S = lyapchol(sys.A,sys.B,sys.E);
        else
            S = lyapchol(sys.A,sys.B);
        end
        S = S';
        if nargout == 2 && Opts.infoLyap == 1 
            data.Info_S = 'No data available for Hammarling method';
        end

        if (nargout > 1 && Opts.infoLyap == 1) || (nargout == 2 && Opts.infoLyap == 0) || nargout == 3 
            if sys.isDescriptor
                R = lyapchol(sys.A', sys.C',sys.E');
            else
                R = lyapchol(sys.A',sys.C');
            end
            R = R';
            if nargout == 3
                data.Info_R = 'No data available for Hammarling method';
            end
        end 
        
    case 'crksm'
        if strcmp(Opts.subspace,'block')
            % get shifts
%             Opts.adiShiftsMethod = 'heur'; % Opts.strategy = 'ADI' / Opts.initShifts??
            [s0_inp,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
            % call CRKSM for S
            [~,V,W,S,dataS] = crksm(sys,s0_inp,Opts);
            if nargout == 2 && Opts.infoLyap == 1 % usage: [S,data] = lyapchol(sys) 
                data.Info_S = dataS;
                data.Info_S.Basis_V = V;
                data.Info_S.Basis_W = W;
            elseif (nargout>1 && Opts.infoLyap == 1) || (nargout == 2 && Opts.infoLyap == 0) || nargout == 3
                % usage: [S,data] = lyapchol(sys), [S,R] = lyapchol(sys), [S,R,data] = lyapchol(sys)
                % call CRKSM for S and R
                if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                    R = S;
                    data.Info_R = dataS; 
                    data.Info_R.Basis_V = V;
                    data.Info_R.Basis_W = W;
                else
                    [~,~,~,R,~] = crksm(sys,[],s0_out,Opts);
                end
            elseif nargout == 3 % usage: [S,R,data] = lyapchol(sys) 
                % call CRKSM for S and R
                if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                    R = S;
                    data.Info_R = dataS; 
                    data.Info_R.Basis_V = V;
                    data.Info_R.Basis_W = W;
                else
                    [~,V,W,R,dataR] = crksm(sys,[],s0_out,Opts);
                    data.Info_R = dataR; 
                    data.Info_R.Basis_V = V;
                    data.Info_R.Basis_W = W;
                end
            end
        else % Opts.subspace='tangential'
            % get shifts and tangential directions
            [s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
            % call CRKSM for S
            [~,V,W,S,dataS] = crksm(sys,s0_inp,Rt,Opts);
            if nargout == 2  && Opts.infoLyap == 1
                data.Info_S = dataS;
                data.Info_S.Basis_V = V;
                data.Info_S.Basis_W = W;
            elseif (nargout>1 && Opts.infoLyap == 1) || (nargout == 2 && Opts.infoLyap == 0) || nargout == 3
                % call CRKSM for S and R
                if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                    R = S;
                    data.Info_R = dataS; 
                    data.Info_R.Basis_V = V;
                    data.Info_R.Basis_W = W;
                else
                    [~,~,~,R,~] = crksm(sys,[],s0_out,[],Lt,Opts);
                end
            elseif nargout == 3
                % call CRKSM for S and R
                if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                    R = S;
                    data.Info_R = dataS;
                    data.Info_R.Basis_V = V;
                    data.Info_R.Basis_W = W;
                else
                    [~,V,W,R,dataR] = crksm(sys,[],s0_out,[],Lt,Opts);
                    data.Info_R = dataR;
                    data.Info_R.Basis_V = V;
                    data.Info_R.Basis_W = W;
                end
            end
        end
        
    otherwise 
        error('sss:lyapchol:invalidMethod','The chosen method for lyapchol is invalid.')
end
% create output
if nargout == 2 && Opts.infoLyap == 1
    varargout{1} = data;
elseif nargout == 2 && Opts.infoLyap == 0
    varargout{1} = R;
elseif nargout == 3
     varargout{1} = R;
     varargout{2} = data;
end
end


