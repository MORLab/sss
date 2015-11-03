function tmax = decayTime(sys)
% Computes the time period in which a sparse LTI system levels off
% ------------------------------------------------------------------
% tmax = decay_time(sys)
% Input:        * sys: an sss-object containing the LTI system
% Output:       * tmax: time after which system has settled down
% ------------------------------------------------------------------
% This file is part of the MORLAB_GUI, a Model Order Reduction and
% System Analysis Toolbox developed at the
% Institute of Automatic Control, Technische Universitaet Muenchen
% For updates and further information please visit www.rt.mw.tum.de
% ------------------------------------------------------------------
% Authors:      Heiko Panzer (heiko@mytum.de), Sylvia Cremer
% Last Change:  28 Oct 2011
% ------------------------------------------------------------------

[res,p]=residue(sys);

% is system stable?
if any(real(p)>0)
    % yes -> tmax=NaN
    tmax=NaN;
    return
end

tmax=0;
for i=1:sys.p
    for j=1:sys.m
        % how much does each pole contribute to energy flow?
        h2=zeros(size(p));
        for k=1:length(p)
            %we need the siso residual for all poles into on vector
            resIJvec = zeros(1,length(p)); 
            for l = 1:length(p), resIJvec(l) = res{l}(i,j);end
            h2(k)=res{k}(i,j)*sum(resIJvec./(-p(k)-p));
        end
        
        [h2_sorted, I] = sort(real(h2));
        % which pole contributes more than 1% of total energy?
        I_dom = I( h2_sorted > 0.01*sum(h2_sorted) );
        if isempty(I_dom)
            % no poles are dominant, use slowest
            I_dom = 1:length(p);
        end
        % use slowest among dominant poles
        [h2_dom, I2] = sort(abs(real(p(I_dom)))); %#ok<ASGLU>

        % when has slowest pole decayed to 1% of its maximum amplitude?
        tmax=max([tmax, log(100)/abs(real(p(I_dom(I2(1)))))]);
    end
end

% store system to caller workspace
if inputname(1)
    sys.decayTime=tmax;
    assignin('caller', inputname(1), sys);
end
