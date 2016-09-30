function [varargout] = last_solveLse(varargin)
% store information about last solveLse call
persistent lastLse

if nargin>0
    lastLse=varargin{1};
end

varargout{1}=lastLse;



end

