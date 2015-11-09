function varargout = disp(sys)
% DISP - Displays information about a sparse state-space model
%
% Syntax:
%       DISP(sys)
%       varargout = DISP(sys)
%
% Description:
%       DISP(sys) displays information about a sparse state-space model
%
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
% Output Arguments:
%       -varargout: character (char) containing the information about the sss-object
%
% Examples:
%       To display some information about the benchmark "build" (SSS, SISO)
%       use:
%
%> load build.mat
%> sys = sss(A,B,C);
%> disp(sys)
%
% See Also: 
%        display
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
% Authors:      Heiko Panzer, Maria Cruz Varona
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  04 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

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

str = [str  char(10), num2str(sys.n) ' states, ' num2str(sys.m) ...
    ' inputs, ' num2str(sys.p) ' outputs'];

if sys.Ts==0
    str = [str  char(10) 'Continuous-time state-space model.'];
else
    str = [str  char(10) 'Sample time: ' num2str(sys.Ts) ' seconds'];
    str = [str  char(10) 'Discrete-time state-space model.'];
end

if nargout>0
    varargout = {str};
else
    str = strrep(str, char(10), [char(10) '  ']);
    disp(['  ' str char(10)]);
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