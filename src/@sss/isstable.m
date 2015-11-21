function [isstable,spectralAbscissa] = isstable(sys)
% ISSTABLE - Check stability of LTI sss system
%
% Syntax:
%       ISSTABLE(sys)
%       isstable = ISSTABLE(sys)
%       [isstable,spectralAbscissa] = ISSTABLE(sys)
%
% Description:
%       This function determines whether the LTI, sss system "sys" is 
%       asymptotically stable. The computations are meant to avoid 
%       operations on full matrices. However, whenever this is not possible,
%       a warning is issued. 
% 
%       The stability check is done by first computing the eigenvalues with 
%       largest real part ('lr'). Afterwards, the spectral abscissa, i.e.
%       the largest occurring real part is calculated. The system is
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
%       The following code checks if the benchmark 'CDplayer' is
%       asymptotically stable:
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> [issd, spectralAbscissa]=isstable(sys);
%
%       Another example, this time using the benchmark 'rail_5177' (DSSS,
%       MIMO):
%
%> load rail_5177.mat
%> sys=sss(A,B,C,[],E);
%> isstable(sys)
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
% with the Chair of Thermofluid Dynamics, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Sylvia Cremer, Alessandro Castagnotto, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%%  Compute the eigenvalue with largest real part
try
    lambda = eigs(sys,1,'lr',struct('v0',sys.b));
catch err
    if strcmp(err.identifier,'MATLAB:eigs:ARPACKroutineErrorMinus14')
        %eigs did not converge: lower the tolerance
        try
            lambda=eigs(sys,1,'lr',struct('tol',1e-4','v0',sys.b));
        catch
            warning('eigs(..,''lr'') failed to compute the spectral abscissa. Trying with eig. This might take a while...');
            lambda = eig(sys);
            lambda = lambda(~isinf(lambda)); %get only finite eigenvalues
        end
    else
        warning('eigs(..,''lr'') failed to compute the spectral abscissa. Trying with eig. This might take a while...');
        lambda = eig(sys);
        lambda = lambda(~isinf(lambda)); %get only finite eigenvalues
    end
end
spectralAbscissa = max(real(lambda));

%%  Check wether the spectral abscissa is strictly less than zero
if  spectralAbscissa < 0
    if nargout<1, fprintf('The system is asymptotically stable\n');
    else isstable = 1; end
        
    elseif spectralAbscissa == 0
        warning('The system has eigenvalues on the imaginary axis. It might be unstable.'); 
        isstable = NaN;
    else
        if nargout<1, warning('The system is unstable.'); else isstable=0;end
end
end