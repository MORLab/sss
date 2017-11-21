function benchmarksCheck(varargin)
% BENCHMARKSCHECK - Verify standard benchmarks 
% 
% Syntax:
%       BENCHMARKSCHECK
%		BENCHMARKSCHECK('sourcePath','/source/path/for/benchmarks/')
%		BENCHMARKSCHECK('destinationPath','/destination/path/for/benchmarks/')
%		BENCHMARKSCHECK('benchmarksList',cellArrayOfBenchmarks)
%
% Description:
%       BENCHMARKSCHECK verifies that the destination path ('destinationPath')
%       contains all the benchmarks specified in a cell array ('benchmarksList')
%       with the filenames of the benchmarks to be verified. In case a benchmark
%       isn't found at the destination path, it is then copied from the source 
%       path.
%
%       When no input arguments are given to BENCHMARKSCHECK, it verifies that a
%       a standard list of benchmarks is available at sss/benchmarks. Any missing
%       benchmark is downloaded from the webdisk.ads.mwn.de server.
%
% Input Arguments:
%       *Parameter Input Arguments:*
%       -'sourcePath':      Source path in case of missing benchmarks.;
%                           'fullfile' or 'filesep' is recommended for path generation.;
%                           Default is webdisk.ads.mwn.de server.
%       -'destinationPath': Path that will be verified for missing benchmarks.;
%                           Default is the sss/benchmarks folder.
%       -'benchmarksList':  Cell array of strings with file names of benchmarks to be verified.;
%                           Example file name: building.mat.;
%                           Default cell array at the beginning of code.
%
%//Note: user defined paths for source and destination can only be folders, not webservers.
%
% Examples:
%       This code verifies if the standard benchmarks are to be found
%       in '~/test/benchmarks' (Unix path example). Missing benchmarks are
%       then downloaded from the webdisk.ads.mwn.de server:
%
%> benchmarksCheck('destinationPath',fullfile('~','test','benchmarks'))
%
% See Also: 
%       loadBenchmarks
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
% Authors:      Rodrigo Mancilla, Mulham Soudan
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?sss">www.rt.mw.tum.de/?sss</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  02 Oct 2017
% Copyright (c) 2015-2017 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------


%% List of standard benchmarks that should always be in the Benchmarks folder
defaultBenchmarksList = {...
	'CDplayer.mat';...
    'LF10.mat';...
    'SpiralInductorPeec.mat';...
    'beam.mat';...
    'building.mat';...
    'eady.mat';...
    'fom.mat';...
    'gyro.mat';...
    'heat-cont.mat';...
    'iss.mat';...
    'bips98_606.mat';...
    'rail_1357.mat';...
    'rail_5177.mat';...
    'random.mat';...
	};

%% Path of local benchmarks folder
[testFolderPath,~,~] = fileparts(mfilename('fullpath'));
sssFolderPath = fileparts(testFolderPath);
benchmarksFolderPath = [sssFolderPath filesep 'benchmarks' filesep];

%% Path to the Web benchmarks folder
% nasFolderPath = [filesep filesep 'nas.ads.mwn.de' filesep 'tumw' filesep 'rtk' filesep 'MORLab' filesep 'Toolbox' filesep '05_benchmark' filesep'];

webSource = 'https://webdisk.ads.mwn.de/Handlers/AnonymousDownload.ashx?folder=6104a79c&path=MORLab%5csssMOR%5cbenchmarks%5c';

%% Parsing
minInputs = 0;
maxInputs = 6;
narginchk(minInputs,maxInputs)

p = inputParser;

addParameter(p,'sourcePath',webSource,@checkPath); % sourcePath is a web address by default. A non default sourcePath can only be a folder address (no web address)
addParameter(p,'destinationPath',benchmarksFolderPath,@checkPath);
addParameter(p,'benchmarksList',defaultBenchmarksList,@iscellstr);

p.FunctionName = 'benchmarks_check';

parse(p,varargin{:});

%% Shorten variable names
sourcePath = p.Results.sourcePath;
destinationPath = p.Results.destinationPath;
benchmarksList = p.Results.benchmarksList;
isWeb = any(strcmp(p.UsingDefaults,'sourcePath'));% Determine if source path is web address (default) or folder address

%% Check for filesep at the end of folder addresses
if destinationPath(end) ~= filesep
    destinationPath = [destinationPath filesep];
end

if ~isWeb % if sourcePath is a folder address (not the default web address)
    if sourcePath(end) ~= filesep
        sourcePath = [sourcePath filesep];
    end
end

%% Determine which mat-files are already in the local folder
filesInFolder = what(destinationPath); % struct with files of all types
matFilesInFolder = filesInFolder.mat; % cell array with names of mat-files

%% Check if all the standard benchmarks are in the local folder
nBenchmarks = length(benchmarksList);

fprintf('\nBenchmarks source:\t%s\nBenchmarks destination:\t%s\n',sourcePath,destinationPath);

if isempty(matFilesInFolder) % if there are no mat-files in the local folder, then copy all benchmarks
    for iBenchmark = 1:nBenchmarks
        transferfiles(sourcePath,destinationPath,isWeb,benchmarksList,iBenchmark)
    end
else
    for iBenchmark = 1:nBenchmarks
        if ~any(strcmp(benchmarksList{iBenchmark}, matFilesInFolder))
            transferfiles(sourcePath,destinationPath,isWeb,benchmarksList,iBenchmark)
        end
    end
end

fprintf('\n\t----> Check complete!');
% The following location contains all standard benchmarks:\n');
% fprintf('\t----> %s\n\n', destinationPath);
% fprintf('-------------------------------------------------------------------------\n');
fprintf('\n-------------------------------------------------------------------------\n\n');

end

function transferfiles(sourcePath,destinationPath,isWeb,benchmarksList,iBenchmark)

% options = weboptions(Name,Value);

if isWeb
    fprintf('Benchmark %s not found. Download in progress...\n',benchmarksList{iBenchmark});
    if verLessThan('MATLAB' , '8.4')
        fprintf('MATLAB R2014a and earlier!\n\n');
        try
            urlwrite([sourcePath char(benchmarksList(iBenchmark))],[destinationPath char(benchmarksList(iBenchmark))]);
        catch
            fprintf('\n\n');
            error('Download was unsuccesful. Please check your internet connection!')
        end
    else
        fprintf('MATLAB R2014b and newer!\n\n');
        try
            currBenchmark = fullfile(destinationPath,benchmarksList{iBenchmark});
            websave(currBenchmark,[sourcePath benchmarksList{iBenchmark}]);
        catch
            fprintf('\n\n');
            error('Download was unsuccesful. Please check your internet connection!')
        end
    end     
    fprintf('%s succesfully downloaded!\n\n',benchmarksList{iBenchmark});
    
else
    fprintf('Benchmark %s not found. Copy in progress...\n',benchmarksList{iBenchmark});
    
    try
        copyfile([sourcePath benchmarksList{iBenchmark}] , destinationPath);
    catch
        fprintf('\n\n');
        error('Copy was unsuccesful. Please check if your source path is correct and accessible!')
    end
    
    fprintf('%s succesfully copied!\n\n',benchmarksList{iBenchmark});
    
end

end

function TF = checkPath(path)

TF = any(exist(path,'dir'));

if ~TF
    error('Couldn''t access/find path. Make sure the path is correct and reachable (e.g. with your file explorer).')
end

end

