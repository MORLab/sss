function [E,A,B,C]=getrail(k);

%% function [E,A,B,C]=getrail(k) reads in data from rail example
% s.a. Oberwolfach Model Reduction Benchmark Collection
%
% Input:
%   k       number of instance 
%           k=1 Data_Rail/rail_1357
%           k=2 Data_Rail/rail_5177 
%           k=3 Data_Rail/rail_20209
%           k=4 Data_Rail/rail_79841
% Output:
%   E       Sparse Matrix
%   A       Sparse Matrix
%   B       Dense Matrix
%   C       Dense Matrix

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

%% set path
switch k
    case 1
        example = 'rail_1357';
    case 2
        example = 'rail_5177';
    case 3  
        example = 'rail_20209';
    case 4
        example = 'rail_79841';
    otherwise
        error('MESS:error_arguments','k has to be 1, 2, 3 or 4\n');
end

%% check path
path = strcat('Data_Rail/',example);
if(~exist(path,'dir'))
   error('MESS:file_io','%s not found \n',path); 
end

%% read matrices
%example file DEMOS/Data_Rail/rail_1357/rail_1357_c60.A
A= mmread(strcat(path,'/',example,'_c60.A'));
B= mmread(strcat(path,'/',example,'_c60.B'));
C= mmread(strcat(path,'/',example,'_c60.C'));
E= mmread(strcat(path,'/',example,'_c60.E'));

%% convert B and C to dense matrices
B = full(B);
C = full(C);





