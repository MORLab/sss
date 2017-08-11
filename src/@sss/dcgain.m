function k = dcgain(sys)
% DCGAIN - Computes the dcgain of a LTI-system
%
% Syntax:
%   k = DCGAIN(sys)
%
% Description:
%       dcgain(sys) computes the dcgain defined as transfer function value
%       at frequency s=0.
%
% Input Arguments:
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
% 
% Output Arguments:
%       -k: dcgain of the system
%
% Examples:
%       The following code computes the dcgain of the benchmark 'CDplayer'
%       (SSS, MIMO):
%
%> load CDplayer.mat
%> sys=sss(A,B,C);
%> k=dcgain(sys)
%
% See Also:
%       freqresp
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
% Authors:      Alessandro Castagnotto, Lisa Jeschek
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  20 May 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

k=freqresp(sys,0);
