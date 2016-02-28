function Opts = parseOpts(Opts,Def)
% PARSEOPTS - Creates an optional struct after parsing
%
% Syntax:
%       Opts = PARSEOPTS(Opts,Def)
%
% Description:
%       TODO
%
% Input Arguments:
%       *Required Input Arguments:*
%           -Opts: A struct containing optional input arguments
%           -Def: default values for the fieldnames defined in the struct Opts
%
% Output Arguments:
%       -Opts: struct containing the optional arguments after the parsing
%
% Examples:
%       TODO
%
% See Also:
%       TODO
%
% References:
%       TODO
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch sss">sss</a>, a Sparse State-Space and System Analysis 
% Toolbox developed at the Chair of Automatic Control in collaboration
% with the Professur fuer Thermofluiddynamik, Technische Universitaet Muenchen. 
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch sss">sss</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Alessandro Castagnotto
% Email:        <a href="mailto:sss@rt.mw.tum.de">sss@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  11 Nov 2015
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

fnames = fieldnames(Def); %target field names

for k = 1:length(fnames)
    if isstruct(Def.(fnames{k})) %nested structure
        % call parseOpts with the nested structure!
          % make sure the nested structure is defined 
            if ~isfield(Opts,fnames{k})
                Opts.(fnames{k}) = struct;
            end
          % call the function again
            Opts.(fnames{k}) = parseOpts(Opts.(fnames{k}),Def.(fnames{k}));
    else
        if ~isfield(Opts,(fnames{k})) || isempty(Opts.(fnames{k}))
            Opts.(fnames{k}) = Def.(fnames{k});
        end
    end
end