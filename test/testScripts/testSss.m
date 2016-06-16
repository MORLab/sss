function result = testSss(Opts)
% TEST - testing of sss functionality
%
% Description:
%   This file provides test functions for the sss toolbox. All tests
%   contained in the folder 'testScripts' can be executed or the tests can
%   be run seperately by choosing the test-suites.
%
%   Required systems: building, beam, fom, eady, random, LF10 (IMTEK), 
%                     SpiralInductorPeec (IMTEK), rail_1357 (IMTEK)
%------------------------------------------------------------------
% This file is part of sssMOR, a Sparse State Space, Model Order
% Reduction and System Analysis Toolbox developed at the Institute 
% of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit www.rt.mw.tum.de
% For any suggestions, submission and/or bug reports, mail us at
%                   -> sssMOR@rt.mw.tum.de <-
%------------------------------------------------------------------
% Authors:      Lisa Jeschek, Jorge Luiz Moreira Silva
% Last Change:  11 Feb 2016
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
% ------------------------------------------------------------------

import matlab.unittest.TestSuite;

%% Choose benchmarks
% Default benchmarks
Def.cond = 'good'; % condition of benchmarks: 'good','bad','all'
                   % 'bad': LF10, beam, random, SpiralInductorPeec
                   % 'good': all benchmarks that are not 'bad'
Def.minSize = 0; % test benchmarks with sys.n >= minSize
Def.maxSize = 400; % test benchmarks with sys.n <= minSize
Def.number = 3; % choose maximum number of tested benchmarks
Def.loadBench = 1; %loadBenchmarks (for testAll)

% create the options structure
if ~exist('Opts','var') || isempty(Opts)
    Opts = Def;
else
    Opts = parseOpts(Opts,Def);
end

% load benchmarks
if Opts.loadBench==1
    testPath=loadBenchmarks(Opts);
else
    p = mfilename('fullpath'); k = strfind(p, fullfile(filesep,'test')); 
	testPath = p(1:k(end)-1);
end

%% Test specific unittest-files
% Available tests:
suite1=TestSuite.fromFile(fullfile(testPath,'testAppend.m'));
suite2=TestSuite.fromFile(fullfile(testPath,'testBode.m'));
suite3=TestSuite.fromFile(fullfile(testPath,'testBodemag.m'));
suite4=TestSuite.fromFile(fullfile(testPath,'testBodeplot.m'));
suite5=TestSuite.fromFile(fullfile(testPath,'testConnect.m'));
suite6=TestSuite.fromFile(fullfile(testPath,'testConnectSss.m'));
suite7=TestSuite.fromFile(fullfile(testPath,'testDcGain.m'));
suite8=TestSuite.fromFile(fullfile(testPath,'testDecayTime.m'));
suite9=TestSuite.fromFile(fullfile(testPath,'testDemoSss.m'));
suite10=TestSuite.fromFile(fullfile(testPath,'testDiag.m'));
suite11=TestSuite.fromFile(fullfile(testPath,'testEig.m'));
suite12=TestSuite.fromFile(fullfile(testPath,'testEigs.m'));
suite13=TestSuite.fromFile(fullfile(testPath,'testFreqresp.m'));
suite14=TestSuite.fromFile(fullfile(testPath,'testImpulse.m'));
suite15=TestSuite.fromFile(fullfile(testPath,'testIssd.m'));
suite16=TestSuite.fromFile(fullfile(testPath,'testIsstable.m'));
suite17=TestSuite.fromFile(fullfile(testPath,'testMtimes.m'));
suite18=TestSuite.fromFile(fullfile(testPath,'testNorm.m'));
suite19=TestSuite.fromFile(fullfile(testPath,'testPoles.m'));
suite20=TestSuite.fromFile(fullfile(testPath,'testPlus.m'));
suite21=TestSuite.fromFile(fullfile(testPath,'testPzmap.m'));
suite22=TestSuite.fromFile(fullfile(testPath,'testResidue.m'));
suite23=TestSuite.fromFile(fullfile(testPath,'testSigma.m'));
% suite24=TestSuite.fromFile(fullfile(testPath,'testSim.m'));
suite25=TestSuite.fromFile(fullfile(testPath,'testSs.m'));
suite26=TestSuite.fromFile(fullfile(testPath,'testStep.m'));
suite27=TestSuite.fromFile(fullfile(testPath,'testZeros.m'));
suite28=TestSuite.fromFile(fullfile(testPath,'testZpk.m'));

suiteSss=[suite1,suite2,suite3,suite4,suite5,suite6,suite7,suite8,suite9,suite10,...
suite11,suite12,suite13,suite14,suite15,suite16,suite17,suite18,suite19,...
suite20, suite21, suite22, suite23, suite25, suite26, suite27, suite28];


%% Run and show results
result = run(suiteSss);
disp(result);
if Opts.loadBench==1
    delete(fullfile(testPath,'benchmarksSysCell.mat'));
end
