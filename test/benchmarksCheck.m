function benchmarksCheck(varargin)



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
    'peec.mat';...
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

webSource = 'https://webdisk.ads.mwn.de/Handlers/AnonymousDownload.ashx?folder=7d2d07de&path=MORLab%5csssMOR%5cbenchmarks%5c';

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

fprintf(2,'\nBenchmarks source:\t%s\nBenchmarks destination:\t%s\n\n',sourcePath,destinationPath);

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

fprintf(2,'-------------------------------------------------------------------------\n');
fprintf(2,'-------------------------------------------------------------------------\n\n');
fprintf(2,'Check complete! The following location contains all standard benchmarks:\n\n%s\n\n', destinationPath);


end

function transferfiles(sourcePath,destinationPath,isWeb,benchmarksList,iBenchmark)

% options = weboptions(Name,Value);

if isWeb
    fprintf(2,'Benchmark %s not found. Download in progress...\n',benchmarksList{iBenchmark});
    
    try
        downloadFile = websave(benchmarksList{iBenchmark},[sourcePath benchmarksList{iBenchmark}]);
    catch
        fprintf('\n\n');
        error('Download was unsuccesful. Please check your internet connection!')
    end
    
    movefile(downloadFile,destinationPath)
    
    fprintf(2,'%s succesfully downloaded!\n\n',benchmarksList{iBenchmark});
    
else
    fprintf(2,'Benchmark %s not found. Copy in progress...\n',benchmarksList{iBenchmark});
    
    try
        copyfile([sourcePath benchmarksList{iBenchmark}] , destinationPath);
    catch
        fprintf('\n\n');
        error('Copy was unsuccesful. Please check if your source path is correct and accessible!')
    end
    
    fprintf(2,'%s succesfully copied!\n\n',benchmarksList{iBenchmark});
end

end

function TF = checkPath(path)

TF = any(exist(path,'dir'));

if ~TF
    error('Couldn''t access/find path. Make sure the path is correct and reachable (e.g. with your file explorer).')
end

end

