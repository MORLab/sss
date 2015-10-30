function  [varargout] = bode(varargin)
% Plots the bode diagram of an LTI system
% ------------------------------------------------------------------
% [mag, phase, omega] = bode(sys, omega, in, out, options)
% Inputs:       * sys: an sss-object containing the LTI system
%    [optional] * vector of imaginary frequencies to plot at
%               * plot options. see <a href="matlab:help plot">PLOT</a>
% Outputs:      * vector of complex frequency response values
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Stefan Jaensch,
%               Sylvia Cremer, Rudy Eid, Alessandro Castagnotto,
%               Lisa Jeschek
% Last Change:  14 Oct 2015
% ------------------------------------------------------------------


% --------------- EVALUATE OPTIONS ---------------
omegaIndex = cellfun(@isfloat,varargin);
omega = [];
if ~isempty(omegaIndex) && nnz(omegaIndex)
    omega = varargin{omegaIndex};
    varargin{omegaIndex}=[];
end

for i = 1:length(varargin)
    if isa(varargin{i},'sss')
        nSys=i;
    end
end

for i = 1:length(varargin)
    if isa(varargin{i},'sss')
        if ~isempty(omega)
            m = freqresp(varargin{i},omega);
        else
            [m, omega] = freqresp(varargin{i});
        end
        
        if isempty(varargin{i}.Name)
            varargin{i}.Name = inputname(i);
        end
        
        % store system in caller workspace
        if inputname(1)
            assignin('caller', inputname(i), varargin{i});
        end
        
        % create frequency response data model
        varargin{i} = frd(m,omega,varargin{i}.Ts,...
            'InputName',varargin{i}.InputName,'OutputName',varargin{i}.OutputName,...
            'Name',varargin{i}.Name);
        
        % clear omega if input consists of more than one system without omega
        if ~nnz(omegaIndex) && nSys>1
            omega=[];
        end
    end
end

% plot
if nargout== 1
    varargout{1} =  varargin{1};
elseif nargout
    [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = bode(varargin{:});
else
    bode(varargin{:});
end
end