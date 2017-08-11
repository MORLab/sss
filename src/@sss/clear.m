function sys = clear(sys)
% CLEAR - Deletes all state-space dimension related properties
%
% Syntax:
%       sys = CLEAR(sys)
%
% Description:
%       sys = clear(sys) deletes the system matrices and labels of the
%       sparse state-space (sss)-object sys.
%
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
% Output Arguments:
%       -sys: the same system with cleared system matrices and labels
%
% Examples:
%> load building.mat;
%> sys=sss(A,B,C);
%> whos('sys'); %shows the size of sys
%> sys=clear(sys);
%> whos('sys') % shows the size of sys with cleared matrices and labels
%
% See Also:
%       clc
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
% Authors:      Thomas Emmert (emmert@tfd.mw.tum.de)
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

sys.A=[];sys.B=[];sys.C=[];sys.D=[];sys.E=[];
sys.InputName=[];
sys.OutputName=[];
sys.StateName=[];
sys.InputGroup=struct();
sys.OutputGroup=struct();
sys.StateGroup=struct();

end