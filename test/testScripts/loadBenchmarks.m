function [testpath] = loadBenchmarks(Opts)
%% testpath
clear;
p = mfilename('fullpath'); k = strfind(p, fullfile(filesep,'test')); 
testpath = [p(1:k(end)-1), fullfile(filesep,'testScripts')];


% Default benchmarks (if testScript is not run by test.m)
Def.cond = 'good'; % condition of benchmarks: 'good','bad','all'
                   % 'bad': LF10, beam, random, SpiralInductorPeec
                   % 'good': all benchmarks that are not 'bad'
Def.minSize = 0; % test benchmarks with sys.n >= minSize
Def.maxSize = 400; % test benchmarks with sys.n <= minSize
Def.number = 5; % choose maximum number of tested benchmarks
Def.dae = 'all'; % 'all', 'withoutDae', 'onlyDae'

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

%% Load benchmarks
%the directory "benchmark" is in sss
p = mfilename('fullpath'); k = strfind(p, fullfile('test',filesep)); 
pathBenchmarks = [p(1:k-1),'benchmarks'];

badBenchmarks = {'LF10.mat','beam.mat','random.mat',...
    'SpiralInductorPeec.mat'};  

% check if benchmarks are in the local benchmarks folder
benchmarksCheck;

% load files
filesAll=dir(fullfile(pathBenchmarks));
files=cell(length(filesAll),1);
nFiles=1; % count of .mat files
for i=1:length(filesAll)
    [~,~,ext]=fileparts(filesAll(i).name);
    if strcmp(ext,'.mat')
        files{nFiles}=filesAll(i).name;
        nFiles=nFiles+1;
    end
end
files(nFiles:end)=[];

benchmarksSysCell=cell(1,Opts.number);
nLoaded=1; %count of loaded benchmarks
disp('Loaded systems:');

warning('off');
for i=1:length(files)
    if nLoaded<Opts.number+1
        sys = sss(files{i});
        if (size(sys.A,1)<=Opts.maxSize && size(sys.A,1)>=Opts.minSize &&...
            ((strcmp(Opts.dae,'all') ||...
             (strcmp(Opts.dae,'withoutDae') && ~sys.isDae)||...
             (strcmp(Opts.dae,'onlyDae') && sys.isDae))))
            switch(Opts.cond)
                case 'good'
                     if ~any(strcmp(files{i},badBenchmarks))
                        if strcmp(files{i},'iss.mat')
                            benchmarksSysCell{nLoaded}=sys([1 3],[1 2 3]);
                        else
                            benchmarksSysCell{nLoaded}=sys;
                        end
                        nLoaded=nLoaded+1;
                        disp(files{i});
                     end
                case 'bad'
                     if any(strcmp(files{i},badBenchmarks))
                        benchmarksSysCell{nLoaded}=sys;
                        nLoaded=nLoaded+1;
                        disp(files{i});
                     end
                case 'all'
                      benchmarksSysCell{nLoaded}=sys;
                      nLoaded=nLoaded+1;
                      disp(files{i});
                otherwise 
                      error('Benchmark option is wrong.');
            end 
        end
    end
end
benchmarksSysCell(nLoaded:end)=[];
warning('on');


% save loaded systems
save(fullfile(fullfile(testpath),'benchmarksSysCell.mat'));
end