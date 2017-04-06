classdef sssTest < matlab.unittest.TestCase
    
    properties 
        sysCell
        deleteBenchmarks
        testPath
    end
    
    methods
        function testCase=sssTest()   
            warning('off','sss:loadSss:secondOrder');
            
            p = mfilename('fullpath'); k = strfind(p, fullfile(filesep,'test')); 
            testCase.testPath = [p(1:k(end)-1), fullfile(filesep,'testScripts')];
            
            if exist(fullfile(testCase.testPath,'benchmarksSysCell.mat'),'file')
                testCase.deleteBenchmarks=0;
            else
                testCase.testPath=loadBenchmarks;
                testCase.deleteBenchmarks=1;
            end
            
            temp=load(fullfile(testCase.testPath,'benchmarksSysCell.mat'));
            testCase.sysCell=temp.benchmarksSysCell;
            if isempty(testCase.sysCell)
                error('No benchmarks loaded.');
            end
        end
    end
    
    methods(TestClassTeardown)
        function deleteBench(testCase)
            warning('on','sss:loadSss:secondOrder');
            if testCase.deleteBenchmarks
                delete(fullfile(testCase.testPath,'benchmarksSysCell.mat'));
            end
        end
    end
end