function  [varargout] = bode(varargin)
% BODE Plots the bode diagram of an LTI system
%
% Syntax:
%   BODE(sys)
%   BODE(sys,omega)
%   BODE(sys1, sys2, ..., omega)
%   BODE(sys1,'-r',sys2,'--k',w);
%   [mag, phase, omega] = BODE(sys)
%   frdData = BODE(sys,omega,'frd')
%
% Description:
%       This function computes the bode plot of one or several LTI systems
%       given as sss objects. If no ouput is specified, then a plot is
%       generated.
%
%       If the frequency range is not specified, the function will
%       determine it. It is also possible to pass several systems or
%       plotting options, just like in MATLAB's built-in version.
%
%       It the function is called with only one ouput and the option 'frd'
%       is specified as last input variable, than an frd object is
%       returned.
%
% Input Arguments:
%       -sys: sss-object containing the LTI system
%       -omega: a vector of frequencies
%
% Output Arguments:
%       - mag/phase: magnitude and phase response
%       - omega:     frequencies corresponding to the data
%       - frdData:   a frd object with the frequency response data
%
% See also:
%   freqresp
%
% References:
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sssMOR">sssMOR</a>, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Chair of
% Automatic Control, Technische Universitaet Muenchen. For updates
% and further information please visit <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sssMOR">sssMOR</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Heiko Panzer, Stefan Jaensch,Sylvia Cremer,
%               Rudy Eid, Alessandro Castagnotto, Lisa Jeschek
% Email:        <a href="mailto:sssMOR@rt.mw.tum.de">sssMOR@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/">www.rt.mw.tum.de</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  01 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

omega = [];
omegaIndex = cellfun(@isfloat,varargin);
if ~isempty(omegaIndex) && nnz(omegaIndex)
    omega = varargin{omegaIndex};
    varargin(omegaIndex)=[];
end

for i = 1:length(varargin)
    % Set name to input variable name if not specified
    if isprop(varargin{i},'Name')
        if isempty(varargin{i}.Name) % Cascaded if is necessary && does not work
            varargin{i}.Name = inputname(i);
        end
    end
    
    % Convert sss to frequency response data model
    if isa(varargin{i},'sss')
        if isempty(omega)
            [m, w] = freqresp(varargin{i});
        else
            [m, w] = freqresp(varargin{i},omega);
        end
        
        %  remove frequencies at infinity to create frd object
        k = find(isinf(w));
        w(k) = []; m(:,:,k) = [];
        
        varargin{i} = frd(m,w,varargin{i}.Ts,...
            'InputName',varargin{i}.InputName,'OutputName',varargin{i}.OutputName,...
            'Name',varargin{i}.Name);
    end
end

% Call ss/bode
if nargout
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = bode(varargin{:});
else
    bode(varargin{:});
end
end