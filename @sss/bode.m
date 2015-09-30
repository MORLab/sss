function varargout = bode(varargin)
% Bode plot for sss objects. 
% ------------------------------------------------------------------
% [mag,phase,wout,sdmag,sdphase] = bode(sys1,sys2,...,sysN,omega,PlotStyleN)
% Inputs:       * sys1: must be sss
%               * omega = 2*pi*f must be specified
% Outputs:      * [mag,phase,wout,sdmag,sdphase] see matlab docu
%               * idfrd object if number of outputs is one
% ------------------------------------------------------------------
% Authors:      Stefan Jaensch, Heiko Panzer
%
% see also: freqresp

omegaIndex = cellfun(@isfloat,varargin);
if all(~omegaIndex)
   error('Frequency vector must be specified for sss/bode') 
end
omega = varargin{omegaIndex};

for i = 1:length(varargin)
    if isa(varargin{i},'sss')
        if varargin{i}.Ts == 0
            G = freqresp(varargin{i},1i* omega);
        else
            G = freqresp(varargin{i},exp(1i* omega*varargin{i}.Ts));
        end        
        if isempty(varargin{i}.Name)
             varargin{i}.Name = inputname(i);
        end        
        varargin{i} = frd(G,omega,varargin{i}.Ts,...
            'InputName',varargin{i}.InputName,'OutputName',varargin{i}.OutputName,...
            'Name',varargin{i}.Name);
    elseif isa(varargin{i},'numlti') || isa(varargin{i},'idlti')
        if isempty(varargin{i}.Name)
             varargin{i}.Name = inputname(i);
        end    
    end    
end
   
if nargout== 1
   varargout{1} =  varargin{1};
elseif nargout
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = bode(varargin{:});    
else
    bode(varargin{:});
end

end