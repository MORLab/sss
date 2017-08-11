function tmax = decayTime(sys)
% DECAYTIME - Computes the time period in which a sparse LTI system levels off
%
% Syntax:
%       DECAYTIME(sys)
%       tmax = DECAYTIME(sys)
%
% Description:
%       tmax = DECAYTIME(sys) computes the time tmax in which the sparse LTI
%       system sys levels off. This is done by taking the slowest pole among
%       the dominant ones and then by computing the time in which the
%       slowest pole decays to 1% of its maximum amplitude.
%
%       If sys is not stable (i.e. there exists at least one pole whose
%       real part is >0), then the decay time is set to tmax=NaN and a 
%       warning is displayed.
%
% Input Arguments:
%       -sys: an sss-object containing the LTI system
%
% Output Arguments:
%       -tmax: time after which the system has settled down
%
% Examples:
%       This code computes the decay time of the benchmark 'building':
%
%> load building; 
%> sys = sss(A,B,C);
%> tmax = decayTime(sys)
%
%       You can visualize the meaning of decay time in a step response plot:
%
%> step(sys);
%> hold on; plot(tmax*[1,1],[-8,8]*1e-4,'r');
%
% See Also:
%       residue, step, impulse
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
% Authors:      Heiko Panzer, Sylvia Cremer, Maria Cruz Varona, 
%               Alessandro Castagnotto
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  19 Jan 2016
% Copyright (c) 2015, 2016 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

%% Compute poles/residues of poles with smallest magnitude
[res,p]=residue(sys,struct('nEigs',min([20,sys.n])));

% is system stable?
if any(real(p)>0 & real(p)<1e6) % larger than 0 but, smaller than infinity-threshold
    % no -> tmax=NaN
    tmax=NaN;
    warning('sss:decayTime:UnstableSys','The system is not stable. The decay time is set to tmax=NaN.');
    return
end

tmax=0; temp = cat(3,res{:}); 
for i=1:sys.p
    for j=1:sys.m
        % how much does each pole contribute to energy flow?
        h2=zeros(size(p));
        for k=1:length(p)
            %we need the siso residual for all poles into on vectors
            resIJvec = squeeze(temp(i,j,:)).';
            h2(k)=res{k}(i,j)*sum(resIJvec./(-p(k)-p).');
        end
        
        [h2_sorted, I] = sort(real(h2));
        % which pole contributes more than 1% of total energy?
        I_dom = I( h2_sorted > 0.01*sum(h2_sorted) );
        if isempty(I_dom)
            % no poles are dominant, use slowest
            I_dom = 1:length(p);
        end
        % use slowest among dominant poles
        [h2_dom, I2] = sort(abs(real(p(I_dom))));

        % when has slowest pole decayed to 1% of its maximum amplitude?
        tmax=max([tmax, log(100)/abs(real(p(I_dom(I2(1)))))]);
    end
end
