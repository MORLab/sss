function [max_r,ind] = mess_s(p,set)
%
% Computation of the maximal magnitude of the rational ADI function over
% a discrete subset of the left complex half plane.
%
%   Calling sequence:
%
%     [max_r,ind] = mess_s(p,set)
%
%   Input:
%
%     p        vector of ADI parameters;
%     set      vector representing the discrete set.
%
%   Output:
%
%     max_r    maximal magnitude of the rational ADI function over set;
%     ind      index - maximum is attained for set(ind). 
%

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
%

%   Exact copy from
%   
%   LYAPACK 1.0 (Thilo Penzl, Jan 1999)
%
if ~isnumeric(p)
    error('MESS:error_arguments','p has to be a vector of numeric type')
end
if ~isnumeric(set)
    error('MESS:error_arguments','set has to be a vector of numeric type')
end
max_r = -1;
ind = 0;
  
for i = 1:length(set)
  
  x = set(i);
  
  rr = 1;
  for j = 1:length(p)

    rr = rr*abs(p(j)-x)/abs(p(j)+x);
    
  end  
    
  if rr > max_r
    
    max_r = rr;
    ind = i;
   
  end
  
end  



