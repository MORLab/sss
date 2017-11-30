function [S,varargout] = lyapchol(varargin)
% LYAPCHOL - Solve Lyapunov equations
%
% Syntax:
%       S               = LYAPCHOL(sys)
%       [S,data]        = LYAPCHOL(sys)
%       [S,R]           = LYAPCHOL(sys)
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
%       If the function is called with [S,R]=lyapchol(sys), then the 
%       low-rank factor Y = R*R' of the dual (generalized) Lyapunov equation 
%       A'*Y*E+E'*Y*A+C'*C=0 is solved as well.
%
%       To call this version of lyapchol with matrices A,B,E, make sure to
%       use |sssFunc.lyapchol(...)|.
%
%       If the option 'method' is set to 'adi', then a low-rank approximation 
%       of the Cholesky (like) factor is performed [2,3]. 
%       If the option 'method' is set to 'crksm', then the Lyapunov equations
%       are approximately solved using CRKSM [1].
%       If this option is not specified, then 'method' is set to 'auto'. In
%       this case, ADI/CRKSM is applied to systems with n>500. For systems
%       with n<=500, 'hammarling' is applied.
%
%       To obtain a data-struct data containing information about the evolution
%       of the residual and some more computed things during the ADI or CRKSM,
%       then set the Opts.infoLyap = 1.
%
%       //Note: the definition of the Cholesky factors X = S*S' is
%       different from built-in lyapchol, where X = S'*S. However, our
%       definition is consistent both with standard literature (cp. [4]) and
%       the low-rank approximation if R has fewer columns than rows.
%
% Input Arguments:
%		*Required Input Arguments:*
%       -sys:   an sss-object containing the LTI system
%       -A,B,E: (alternatively) matrices of the Lyapunov equation
%		*Optional Input Arguments:*
%       -Opts:                     a structure containing following options
%           -.method:              select solver for lyapunov equation 
%                                  [{'auto'} / 'adi' / 'hammarling' / 'crksm']
%           -.lse:                 solve linear system of equations (only ADI and CRKSM)
%                                  [{'gauss'} / 'luChol']
%           -.rctol:               tolerance for difference between ADI or CRKSM iterates
%                                  [{'1e-12'} / positive float]
%           -.restol:              tolerance for the residual of the Lyapunov equation
%                                  [{'1e-8'} / positive float]
%           -.maxiter:             maximum number of iterations (only ADI and CRKSM)
%                                  [{150} / positive integer]
%           -.infoLyap:            output data-struct in varargout no/yes?
%                                  [{0} / 1]
%           -.adi.norm             refer to opts.adi.norm in MESS_LRADI for 
%                                  more info  
%                                  [2, {'fro'}]
%           -.adi.shifts.l0:       refer to opts.adi.shifts.l0 in
%                                  MESS_PARA or MESS_LRADI for more info
%                                  [{20} / positive integer]
%           -.adi.shifts.kp:       refer to opts.adi.shifts.kp in
%                                  MESS_PARA or MESS_LRADI for more info
%                                  [{50} / positive integer]
%           -.adi.shifts.km:       refer to opts.adi.shifts.km in
%                                  MESS_PARA or MESS_LRADI for more info
%                                  [{25} / positive integer]
%           -.adi.shifts.method:   refer to opts.adi.shifts.method in
%                                  MESS_PARA or MESS_LRADI for more info
%                                  [{'heur'} / 'wachspress' / 'projection']
%           -.q:                   size of Cholesky factor (only ADI)
%                                  [{'0'} / positive integer]
%           -.forceOrder:          return order q
%                                  [{'false'} / 'true']
%           -.subspace:            build block or tangential Krylov subspace (only CRKSM)
%                                  [{'block'} / 'tangential']
%           -.initShiftsStrategy:  choose initial shift strategy in INITIALIZESHIFTS (only CRKSM)
%                                  [{'ADI'} / 'eigs' / 'ROM' / 'const']
%           -.nShifts:             set number of initial shifts in INITIALIZESHIFTS (only CRKSM)
%                                  [{10} / positive, even integer]
%
% Output Arguments:
%       -S:     Cholesky factor X=S*S' of Lyapunov equation A*X*E'+E*X*A'+B*B'=0
%       -R:     Cholesky factor Y=R*R' of Lyapunov equation A'*Y*E+E'*Y*A+C'*C=0
%
% Examples:
%       Compute the Cholesky factors for both Lyapunov equations
%
%> sys = sss('building');
%> [S,R] = lyapchol(sys); % using 'hammarling', since n<=500
%
%       To compute a single Cholesky factor, use
%
%> sys = sss('fom'); Opts.infoLyap = 1;
%> [Sadi,dataSadi] = lyapchol(sys,Opts); % using 'auto', since n<=500
%> [Radi,dataRadi] = lyapchol(sys',Opts);
%
%       Compute the Cholesky factors for both Lyapunov equations using CRKSM
%
%> clear
%> sys = sss('fom'); Opts.method = 'crksm';
%> [Scrksm,Rcrksm,dataCrksm] = lyapchol(sys,Opts);
%> ScrksmLR = dataCrksm.Info_S.Basis_V*Scrksm'; % compute low-rank factor
%
%       To call this version of lyapchol with matrices use
%
%> sys = sss('building');
%> [A,B,C,D,E] = dssdata(sys);
%> S = sssFunc.lyapchol(A,B,E);
%
% See Also:
%       solveLse, tbr, norm, numerics/lyapchol, mess_lradi, crksm
%
% References:
%       * *[1] Druskin, Simoncini (2011)*, Adaptive Rational Krylov Subspaces
%       for large-scale dynamical systems
%       * *[2] Penzl (2000)*, LYAPACK - A MATLAB Toolbox for Large Lyapunov
%       and Riccati Equations, Model Reduction Problems, and Linear-Quadratic 
%       Optimal Control Problems.
%       * *[3] Saak, Koehler, Benner (2016)*, M-M.E.S.S. - The Matrix 
%       Equation Sparse Solver Library.
%       * *[4] Golub, Van Loan (1996)*, Matrix computations
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
% Last Change:  23 Nov 2017
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
Def.method              = 'auto';           % 'auto', 'adi', 'hammarling', 'crksm'
Def.lse                 = 'gauss';          % only for MESS (see solveLse) and CRKSM
Def.restol              = 1e-8;             % only for MESS and CRKSM
Def.rctol               = 0;                % only for MESS and CRKSM
Def.maxiter             = min([200,sys.n]); % only for MESS and CRKSM
Def.infoLyap            = 0;                % output data-struct in varargout no/yes?

% ADI default execution parameters
Def.adi.norm            = 'fro';            % only for MESS
Def.adi.shifts.l0       = 20;               % only for MESS
Def.adi.shifts.kp       = 50;               % only for MESS
Def.adi.shifts.km       = 26;               % only for MESS
Def.adi.shifts.method   = 'heur';           % only for MESS  
Def.q                   = 0;                % only for MESS
Def.forceOrder          = false;            % only for MESS

% CRKSM default execution parameters
Def.subspace            = 'block';          % only for CRKSM; build block or tangential Krylov subspace, [{'block'} / 'tangential']
Def.initShiftsStrategy  = 'ADI';            % only for CRKSM; choose shift strategy in INITIALIZESHIFTS, [{'ADI'} / 'eigs' / 'ROM' / 'const']
Def.nShifts             = 10;               % only for CRKSM; set number of initial shifts

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
        Opts.method = 'adi';    %set ADI/CRKSM as default in case 'auto'
    else
        Opts.method = 'hammarling';
    end
end

 % check if CRKSM and internal functions are available
if strcmp(Opts.method, 'crksm') && (~exist('crksm.m','file') || ~exist('initializeShifts.m','file') || ~exist('getShifts.m','file')...
   || ~exist('solveLse.m','file')) 
    warning(['One of the following files not found:',...
        'crksm.m, initializeShifts.m, getShifts.m, solveLse.m',...
        'Using ADI method instead of CRKSM.']);
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
        messOpts.adi=struct('shifts',struct('l0',Opts.adi.shifts.l0,'kp',Opts.adi.shifts.kp,...
            'km',Opts.adi.shifts.km,'b0',ones(sys.n,1),'method',Opts.adi.shifts.method,'info',0),...
            'maxiter',Opts.maxiter,'restol',Opts.restol,'rctol',Opts.rctol,...
            'info',0,'norm',Opts.adi.norm);

        oper = operatormanager(lseType);
        messOpts.solveLse.lse=Opts.lse;
        messOpts.solveLse.krylov=0;

        % get adi shifts
        [messOpts.adi.shifts.p]=mess_para(eqn,messOpts,oper);

        % low-rank adi
        [S,Sout]=mess_lradi(eqn,messOpts,oper);

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
        
        if nargout == 2 && Opts.infoLyap == 1 % usage: [S,data] = lyapchol(sys) 
            data.Info_S = Sout;
        elseif (nargout == 2 && Opts.infoLyap == 0) || nargout == 3 % usage: [S,R] = lyapchol(sys), [S,R,data] = lyapchol(sys)
            if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                R=S;
            else
                eqn.type='T';
                [R,Rout]=mess_lradi(eqn,messOpts,oper);
                
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
            if nargout == 3 % usage: [S,R,data] = lyapchol(sys)
                if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                    data.Info_S = Sout;
                    data.Info_R = data.Info_S;
                else
                    data.Info_S = Sout;     data.Info_R = Rout;
                end
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
        
        if nargout == 2 && Opts.infoLyap == 1 % usage: [S,data] = lyapchol(sys)
            data.Info_S = 'No data available for Hammarling method';
        elseif (nargout == 2 && Opts.infoLyap == 0) || nargout == 3 % usage: [S,R] = lyapchol(sys), [S,R,data] = lyapchol(sys)
            if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                R=S';
            else
                if sys.isDescriptor
                    R = lyapchol(sys.A',sys.C',sys.E');
                else
                    R = lyapchol(sys.A',sys.C');
                end
            end
            R = R';
            if nargout == 3 % usage: [S,R,data] = lyapchol(sys)
                data.Info_S = 'No data available for Hammarling method';
                data.Info_R = 'No data available for Hammarling method';
            end
        end
        
    case 'crksm'
        if strcmp(Opts.subspace,'block')
            % get shifts
            [s0_inp,~,s0_out] = initializeShifts(sys,Opts.nShifts,1,Opts);
            % call CRKSM for S
            [~,V_S,W_S,S,dataS] = crksm(sys,s0_inp,Opts);
            
            if nargout == 2 && Opts.infoLyap == 1 % usage: [S,data] = lyapchol(sys) 
                data.Info_S = dataS;
                data.Info_S.Basis_V = V_S;
                data.Info_S.Basis_W = W_S;
            elseif (nargout == 2 && Opts.infoLyap == 0) || nargout == 3 % usage: [S,R] = lyapchol(sys), [S,R,data] = lyapchol(sys)
                % call CRKSM for R
                if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                    R = S;
                else
                    [~,V_R,W_R,R,dataR] = crksm(sys,[],s0_out,Opts);
                end
                if nargout == 3 % usage: [S,R,data] = lyapchol(sys)
                    if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                        data.Info_S = dataS;           
                        data.Info_S.Basis_V = V_S;     
                        data.Info_S.Basis_W = W_S;
                        data.Info_R = data.Info_S;
                    else
                        data.Info_S = dataS;           data.Info_R = dataR;
                        data.Info_S.Basis_V = V_S;     data.Info_R.Basis_V = V_R;
                        data.Info_S.Basis_W = W_S;     data.Info_R.Basis_W = W_R;
                    end
                end
            end
        else % Opts.subspace='tangential'
            % get shifts and tangential directions
            Opts.initShiftsStrategy = 'eigs';
            fprintf('initial shift strategy is automatically set to "eigs" because ADI strategy does not support the tangential case');
            [s0_inp,Rt,s0_out,Lt] = initializeShifts(sys,Opts.nShifts,1,Opts);
            % call CRKSM for S
            [~,V_S,W_S,S,dataS] = crksm(sys,s0_inp,Rt,Opts);
            
            if nargout == 2 && Opts.infoLyap == 1 % usage: [S,data] = lyapchol(sys)
                data.Info_S = dataS;
                data.Info_S.Basis_V = V_S;
                data.Info_S.Basis_W = W_S;
            elseif (nargout == 2 && Opts.infoLyap == 0) || nargout == 3 % usage: [S,R] = lyapchol(sys), [S,R,data] = lyapchol(sys)
                % call CRKSM for R
                if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                    R = S; 
                else
                    [~,V_R,W_R,R,dataR] = crksm(sys,[],s0_out,[],Lt,Opts);
                end
                if nargout == 3 % usage: [S,R,data] = lyapchol(sys)
                    if sys.isSym && ~any(size(sys.B)-size(sys.C')) && norm(full(sys.B-sys.C'))==0
                        data.Info_S = dataS;           
                        data.Info_S.Basis_V = V_S;     
                        data.Info_S.Basis_W = W_S;
                        data.Info_R = data.Info_S;
                    else
                        data.Info_S = dataS;           data.Info_R = dataR;
                        data.Info_S.Basis_V = V_S;     data.Info_R.Basis_V = V_R;
                        data.Info_S.Basis_W = W_S;     data.Info_R.Basis_W = W_R;
                    end
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
