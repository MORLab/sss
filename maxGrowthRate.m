function stable = maxGrowthRate(sys)
% maxGrowthRate function determines the maximum growth rate
% ------------------------------------------------------------------
% This file is part of tax, a code designed to investigate thermoacoustic
% network systems. It is developed by:
% Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen.
% For updates and further information please visit www.tfd.mw.tum.de
% ------------------------------------------------------------------
% stable = stability(sys);
% Input:        * sys: tax object
% Output:       * stable: maximum growth rate. 
%                         Continuous time: maximum real part;
%                         Discrete time: maximum absolute value.
% ------------------------------------------------------------------
% Authors:      Stefan Jaensch (jaensch@tfd.mw.tum.de)
% Last Change:  23 Out 2015
% ------------------------------------------------------------------

if sys.Ts == 0
    stable = eigs(sys,'lr');
else
    stable = eigs(sys,'lm');
end