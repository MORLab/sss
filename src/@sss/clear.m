function sys = clear(sys)
% CLEAR - Deletes all state-space dimension related properties
%
% Syntax:
%       sys = CLEAR(sys)
%
% Description:
%       sys = clear(sys) deletes the system matrices and lables of the
%       sparse state-space (sss)-object sys.
%
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
% Output Arguments:
%       -sys: the same system with cleared system matrices and labels
%
% Examples:
%       TODO
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
% Authors:      Thomas Emmert (emmert@tfd.mw.tum.de)
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
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