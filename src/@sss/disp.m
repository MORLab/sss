function varargout = disp(sys)
%% Display functions
mc = metaclass(sys);
str = [];
if ~isempty(mc.Name) && ~isempty(sys.Name)
    str = [mc.Name ' Model ' sys.Name];
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

if ~isempty(sys.morInfo)
    str = [str char(10) sys.dispMorInfo];
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