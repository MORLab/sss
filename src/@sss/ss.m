function sys_ss = ss(sys_sss)
% SS - Converts sparse LTI system (sss) to Matlab\control\ss
%
% Syntax:
%       sys_ss = ss(sys_sss)
%
% Description:
%       Converts sparse LTI system (sss) to Matlab\control\ss
%
% Input Arguments:
%       -sys_sss: sparse state space (sss)-object
%
% Output Arguments:       
%       -sys_ss:  ss- or dss-object, respectively
%
% Examples:
%       TODO
%
% See Also:
%       TODO
%
% References:
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
% Authors:      Heiko Panzer, Sylvia Cremer, Thomas Emmert
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  05 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if sys_sss.isDescriptor
    sys_ss=dss(full(sys_sss.A),full(sys_sss.B),full(sys_sss.C),full(sys_sss.D),full(sys_sss.E));
else
    sys_ss=ss(full(sys_sss.A),full(sys_sss.B),full(sys_sss.C),full(sys_sss.D));
end
sys_ss.Ts = sys_sss.Ts;

sys_ss.y = sys_sss.y;
sys_ss.u = sys_sss.u;
sys_ss.StateName = sys_sss.StateName;
sys_ss.Name = sys_sss.Name;

for field = fieldnames(sys_sss.OutputGroup)'
    sys_ss.OutputGroup.(char(field)) = sys_sss.OutputGroup.(char(field));
end
for field = fieldnames(sys_sss.InputGroup)'
    sys_ss.InputGroup.(char(field)) = sys_sss.InputGroup.(char(field));
end
for field = fieldnames(sys_sss.StateGroup)'
    sys_ss.UserData.StateGroup.(char(field)) = sys_sss.StateGroup.(char(field));
end
if isstruct(sys_sss.UserData)
    for field = fieldnames(sys_sss.UserData)'
        sys_ss.UserData.(char(field)) = sys_sss.UserData.(char(field));
    end
end

