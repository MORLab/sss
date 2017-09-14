function varargin = getfrd(varargin)
% GETFRD - Get frd-object(s) for frequency-response functions
% 
% Syntax: 
%       varargin = GETFRD(sys)
%       varargin = GETFRD(sys, omega)
%       varargin = GETFRD(sys, LineSpec, omega)
%       varargin = GETFRD(sys1, sys2, ..., sysN)
%       varargin = GETFRD(sys1, sys2, ..., sysN, omega)
%       varargin = GETFRD(sys1,'-r',sys2,'--k');
%       varargin = GETFRD(sys1,'-r',sys2,'--k',omega)
%
% Description:
%       Auxiliary function that is used within |bode|, |bodemag|, etc. to
%       calculate and get the frd-object(s). 
%
% Input Arguments:       
%       *Required Input Arguments:*
%       -sys: an sss-object containing the LTI system
%       *Optional Input Arguments:*
%       -omega: vector of frequencies or cell with {wmin,wmax}
%       -LineSpec: String with line style, marker symbol, and color. See <a href="matlab:doc plot">plot</a>
%
% Output Arguments:      
%       -varargin: cell array with frd-objects and eventually line styles
%
% See Also:
%       bode, bodemag, bodeplot, sigma, freqresp, frd
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
% Authors:      Maria Cruz Varona, Alessandro Castagnotto 
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  14 Sep 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

%% Check for omega
omegaIndex=0;
for i=1:length(varargin)
    if isa(varargin{i}, 'double')|| isa(varargin{i},'cell')
        omegaIndex = i;
        break
    end
end

if omegaIndex>0
	omega=varargin{omegaIndex};
    varargin(omegaIndex)=[];
else
    omega=[];
end

%% Get frd object
for i = 1:length(varargin)
    % Set name to input variable name if not specified
    if isprop(varargin{i},'Name')
        if isempty(varargin{i}.Name) % Cascaded if is necessary && does not work
            varargin{i}.Name = inputname(i);
        end
    end
    
    % Convert sss to frequency response data model
    if isa(varargin{i},'sss') || isa(varargin{i},'ssRed')
        if not(exist('omega','var')) || isempty(omega)
            sysfrd{i}   = frd(varargin{i}); %#ok<*AGROW>
            omegaMin(i) = sysfrd{i}.Frequency(1); %#ok<*AGROW>
            omegaMax(i) = sysfrd{i}.Frequency(end); %#ok<*AGROW>
        else
            varargin{i} = frd(varargin{i},omega);
        end
    else
        if not(exist('omega','var')) || isempty(omega)
            omegaMin(i) = NaN;
            omegaMax(i) = NaN;
        end
    end
end

%% make sure that all plots have the same frequency range if 'omega' is not specified
if not(exist('omega','var')) || isempty(omega)
    wMin = min(omegaMin);
    wMax = max(omegaMax);

    for i = 1:length(varargin)
        if isa(varargin{i},'sss') || isa(varargin{i},'ssRed')
            if wMin < sysfrd{i}.Frequency(1) % adjustment of wMin
                sysfrdwMin = frd(varargin{i},wMin);
                G = cat(3, sysfrdwMin.ResponseData, sysfrd{i}.ResponseData);
                omega = cat(1, wMin, sysfrd{i}.Frequency);
                sysfrd{i} = frd(G,omega,varargin{i}.Ts,...
                'InputName',varargin{i}.InputName,'OutputName',varargin{i}.OutputName,...
                'Name',varargin{i}.Name);
            end

            if wMax > sysfrd{i}.Frequency(end) % adjustment of wMax
                sysfrdwMax = frd(varargin{i},wMax);
                G = cat(3, sysfrd{i}.ResponseData, sysfrdwMax.ResponseData);
                omega = cat(1, sysfrd{i}.Frequency, wMax);
                sysfrd{i} = frd(G,omega,varargin{i}.Ts,...
                'InputName',varargin{i}.InputName,'OutputName',varargin{i}.OutputName,...
                'Name',varargin{i}.Name);
            end
            varargin{i} = sysfrd{i};
        end
    end
end

end
