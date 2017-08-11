function infostr = disp(sys)
% DISP - Displays information about a sparse state-space model
%
% Syntax:
%       DISP(sys)
%       infostr = DISP(sys)
%
% Description:
%       DISP(sys) displays information about a sparse state-space model:
%        
%       # SSS, DSSS or DAE
%       # SISO, SIMO, MISO, MIMO
%       # Number of state, input and output variables
%       # Continous or discrete-time state-space model
% 
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
% Output Arguments:
%       -infostr: string (char) containing the information about the sss-object
%
% Examples:
%       To display some information about the benchmark "iss" (SSS, MIMO)
%       use:
%
%> load iss.mat
%> sys = sss(A,B,C);
%> disp(sys)
%
% See Also: 
%        display
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
% Authors:      Heiko Panzer, Maria Cruz Varona
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

if isempty(sys)
    fprintf(1,'  Empty sparse state-space model.\n\n');
else
    mc = metaclass(sys);
    str = [];
    if ~isempty(mc.Name) && ~isempty(sys.Name)
        str = [mc.Name ' Model ' sys.Name, ' '];
    end

    if sys.isDae;            str = [str '(DAE)'];
    elseif sys.isDescriptor; str = [str '(DSSS)'];
    else                     str = [str '(SSS)'];
    end

    if sys.isSiso;       str = [str '(SISO)'];
    elseif sys.isSimo;   str = [str '(SIMO)'];
    elseif sys.isMiso;   str = [str '(MISO)'];
    elseif sys.isMimo;   str = [str '(MIMO)'];
    end

    str = [str  char(10), num2str(sys.n) ' state variables, ' num2str(sys.m) ...
        ' inputs, ' num2str(sys.p) ' outputs'];

    if sys.Ts==0
        str = [str  char(10) 'Continuous-time state-space model.'];
    else
        str = [str  char(10) 'Sample time: ' num2str(sys.Ts) ' seconds'];
        str = [str  char(10) 'Discrete-time state-space model.'];
    end

    if nargout>0
        infostr = {str};
    else
        str = strrep(str, char(10), [char(10) '  ']);
        disp(['  ' str char(10)]);
    end
end

end

function str = dispMorInfo(sys)
if isempty(sys.morInfo)
    str=[];
else
    str = ['Created by MOR on ' datestr(sys.morInfo.time) '.' char(10) ...
        'Reduction Method: ' sys.morInfo.method  '.' char(10) ...
        'Original system: ' sys.morInfo.orgsys '.'];
end
end