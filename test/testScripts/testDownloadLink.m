classdef testDownloadLink < sssTest
% testDownloadLink - Testing if the download links work
%
% Description:
%   The download link are tested for:
%    + Download of a benchmark model
%
%------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
%------------------------------------------------------------------
% Authors:      Niklas Kochdumper
% Last Change:  19 Jul 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    methods(Test)
        
        function testDownloadLinkBenchmark(testCase)
        %% downlaod benchmark CDplayer and load it to the workspace in 
        %  order to test if the download was sucessful
        
            % download-link 
            webSource = 'https://webdisk.ads.mwn.de/Handlers/AnonymousDownload.ashx?folder=6104a79c&path=MORLab%5csssMOR%5cbenchmarks%5cCDplayer.mat';
            
            % file path to store the file
            pathRoot = fileparts(mfilename('fullpath'));
            path = fullfile(pathRoot,'Temp.mat');
            
            % download the file
            websave(path,webSource);
            
            % load the benchmark
            sys = loadSss('Temp.mat');
            
            % delete the downloaded file 
            delete(path);
            
        end 
    end
end