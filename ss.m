function sys_ss = ss(sys_sss)
% Converts sparse LTI system (sss) to Matlab\control\ss
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% sys_ss = ss(sys_sss);
% Input:        * sys_sss: sparse state space (sss)-object
% Output:       * sys_ss:  ss- or dss-object, respectively
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer,
%               Thomas Emmert (emmert@tfd.mw.tum.de)
% Last Change:  17 Feb 2015
% ------------------------------------------------------------------

if sys_sss.is_dae
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

