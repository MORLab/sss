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
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Sylvia Cremer, Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if sys.n < 100
    %%  For small systems, compute the eigenvalue decomposition directy
    lambda = eig(sys);
else       
    %%  Compute the eigenvalue with largest real part
    try
        lambda = eigs(sys,1,'lr',struct('v0',sys.b));
    catch err
        if strcmp(err.identifier,'MATLAB:eigs:ARPACKroutineErrorMinus14')
            %eigs did not converge: lower the tolerance
            try
                lambda=eigs(sys,1,'lr',struct('tol',1e-4','v0',sys.b));
            catch
                warning('sss:isstable:EigsFailed','eigs(..,''lr'') failed to compute the spectral abscissa. Trying with eig. This might take a while...');
                lambda = eig(sys);
            end
        else
            warning('sss:isstable:EigsFailed','eigs(..,''lr'') failed to compute the spectral abscissa. Trying with eig. This might take a while...');
            lambda = eig(sys);
        end
    end
end
lambda = lambda(~isinf(lambda)); %get only finite eigenvalues
lambda = lambda(abs(real(lambda))<1e6); % infinity-threshold
spectralAbscissa = max(real(lambda));

if isempty(spectralAbscissa)
    warning('The spectral abscissa is empty.');
    isstable=NaN;

%%  Check wether the spectral abscissa is strictly less than zero
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