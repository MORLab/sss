function [isstable,spectralAbscissa] = isstable(sys)
% ISSTABLE - Check stability of sparse LTI system
%
% Syntax:
%       ISSTABLE(sys)
%       isstable = ISSTABLE(sys)
%       [isstable,spectralAbscissa] = ISSTABLE(sys)
%
% Description:
%       This function determines whether the sparse LTI system |sys| is 
%       asymptotically stable. The computations are meant to avoid 
%       operations on full matrices. However, whenever this is not possible,
%       a warning is issued. 
% 
%       The stability check is done by first computing the eigenvalues with 
%       largest real part ('lr'). Afterwards, the spectral abscissa, i.e.
%       the largest occurring real part, is calculated. The system is
%       asymptotically stable, if the spectral abscissa is strictly less 
%       than zero.
% 
%       If no output is defined, then the result is printed on the screen.
%       Depending on the number of ouputs defined the function can return
%       isstable and spectralAbscissa.
% 
%       NaN is returned either, when the computation was not possible, or 
%       when the spectral abscissa is zero. In the latter case, the system
%       might be stable (in the sense of Lyapunov) or unstable if the
%       multiplicity of the eigenvalues at the origin is greater than one.
%
%       //Note: if the computation of the spectral abscissa with eigs
%       fails, then this function uses eig to determine the eigenvalue with
%       largest real part. This is possible only for mid-sized problems.
%
% Input Arguments:
%       -sys: sss-object containing the LTI system
%
% Output Arguments:
%       -isstable: a boolean value (1=true, 0=false, NaN=unknown)
%       -spectralAbscissa: i.e. the largest real part of the eigenvalues.
%
% Examples:
%       Following example determines the stability of the 'CDplayer' model
%       (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> isstable(sys)
%
%      To obtain the spectral abscissa the function can be called with two
%      outputs
%
%> load building.mat
%> sys=sss(A,B,C);
%> [isstab, spectralAbscissa]=isstable(sys)
%
%
% See Also:
%       issd, eigs, chol, pzmap
%
% References:
%       * *[1] Panzer (2014)*, Model order reduction by Krylov subspace methods
%       with Global Error Bounds and Automatic Choice of Parameters
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
% Authors:      Sylvia Cremer, Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  26 Oct 2018
% Copyright (c) 2015-2018 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if sys.n < 100
    %%  For small systems, compute the eigenvalue decomposition directy
    lambda = eig(sys);
else       
    %%  Compute the eigenvalue with largest real part
    if issymmetric(sys.A) && issymmetric(sys.E) % A = A' and E = E'
        eigsOpt = 'la';
        isSymPosDef = ispd(sys.e); % E > 0
    else
        eigsOpt = 'lr'; 
        isSymPosDef = 0;
    end
    [~,lambda,flag] = eigs(sys,1,eigsOpt,struct('spdB',isSymPosDef,'v0',sum(sys.e,2)));
    
    if flag
        %eigs did not converge: lower the tolerance
        [~,lambda,flag]=eigs(sys,1,'lr',struct('tol',1e-4,'spdB',isSymPosDef,'v0',sum(sys.e,2)));
        
        if flag
            %eigs did not converge: try with eig
            warning('sss:isstable:EigsFailed','eigs(..,''lr'') failed to compute the spectral abscissa. Trying with eig. This might take a while...');
            lambda = eig(sys);
        end
    end
end
lambda = lambda(~isinf(lambda)); %get only finite eigenvalues
lambda = lambda(abs(real(lambda))<1e8); % infinity-threshold
spectralAbscissa = max(real(lambda));

if isempty(spectralAbscissa)
    warning('The spectral abscissa is empty.');
    isstable=NaN;

%%  Check whether the spectral abscissa is strictly less than zero
elseif  spectralAbscissa < 0
    if nargout<1, fprintf('The system is asymptotically stable.\n');
    else isstable = 1; end
        
    elseif spectralAbscissa == 0
        warning('The system has eigenvalues on the imaginary axis. It might be unstable.'); 
        isstable = NaN;
    else
        if nargout<1, fprintf('The system is unstable.\n'); else isstable=0;end
end
end