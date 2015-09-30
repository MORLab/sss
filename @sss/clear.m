function sys = clear(sys)
%clear deletes all state space dimension related properties
% ------------------------------------------------------------------
% [sys] = clear(sys)
% Input:        * sys: an sss-object containing the LTI system
% Output:       * sys: the same system with cleared system matrices and
%                      labels
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Thomas Emmert (emmert@tfd.tum.de)
% Last Change:  25 Mar 2015
% ------------------------------------------------------------------

sys.A=[];sys.B=[];sys.C=[];sys.D=[];sys.E=[];
sys.InputName=[];
sys.OutputName=[];
sys.StateName=[];
sys.InputGroup=struct();
sys.OutputGroup=struct();
sys.StateGroup=struct();

end