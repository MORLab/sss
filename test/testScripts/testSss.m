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

%%  Change to testScripts folder
testCase.Path = pwd; %original

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

% load benchmarks and change to folder 'testScripts'
if Opts.loadBench==1
    loadBenchmarks(Opts);
end

%% Test specific unittest-files
% Available tests:
suite1=TestSuite.fromFile('testAppend.m');
suite2=TestSuite.fromFile('testBode.m');
suite3=TestSuite.fromFile('testBodemag.m');
suite4=TestSuite.fromFile('testBodeplot.m');
suite5=TestSuite.fromFile('testConnect.m');
suite6=TestSuite.fromFile('testConnectSss.m');
suite7=TestSuite.fromFile('testDecayTime.m');
suite8=TestSuite.fromFile('testDiag.m');
suite9=TestSuite.fromFile('testEig.m');
suite10=TestSuite.fromFile('testEigs.m');
suite11=TestSuite.fromFile('testFreqresp.m');
suite12=TestSuite.fromFile('testImpulse.m');
suite13=TestSuite.fromFile('testIssd.m');
suite14=TestSuite.fromFile('testIsstable.m');
suite15=TestSuite.fromFile('testMtimes.m');
suite16=TestSuite.fromFile('testNorm.m');
suite17=TestSuite.fromFile('testPlus.m');
suite18=TestSuite.fromFile('testPzmap.m');
suite19=TestSuite.fromFile('testResidue.m');
suite20=TestSuite.fromFile('testSigma.m');
% suite21=TestSuite.fromFile('testSim.m');
suite22=TestSuite.fromFile('testSs.m');
suite23=TestSuite.fromFile('testStep.m');

suiteSss=[suite1,suite2,suite3,suite4,suite5,suite6,suite7,suite8,suite9,suite10,...
suite11,suite12,suite13,suite14,suite15,suite16,suite17,suite18,suite19,...
suite20, suite22, suite23];

%% Run and show results
result = run(suiteSss);
disp(result);
if Opts.loadBench==1
    delete('benchmarksSysCell.mat');
end

%% Go back to original folder
cd(testCase.Path);
