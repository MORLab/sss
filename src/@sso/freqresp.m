function [G,w] = freqresp(sys,w)

%   very simple version until the full one works..
G = zeros(sys.p,sys.m,length(w));
for iW = 1:length(w);
    G(:,:,iW) = (sys.Cp + 1i*w(iW)*sys.Cv)*...
                ((1i*w(iW)^2*sys.M+1i*w(iW)*sys.D+sys.K)\sys.B);
end