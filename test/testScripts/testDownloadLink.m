classdef testDownloadLink < sssTest
% testDownloadLink - Testing if the download links work
%
% Description:
%   The download link are tested for:
%    + Download of a benchmark model
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
%------------------------------------------------------------------
% Authors:      Niklas Kochdumper, Mulham Soudan, Maria Cruz Varona
% Last Change:  02 Oct 2017
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
            if verLessThan('MATLAB' , '8.4')
                fprintf('MATLAB R2014a and earlier!\n\n');
                urlwrite(webSource,path);
            else
                fprintf('MATLAB R2014b and newer!\n\n');
                websave(path,webSource);
            end
            
            % load the benchmark
            sys = sss('Temp.mat');
            
            % delete the downloaded file 
            delete(path);
            
        end 
    end
end