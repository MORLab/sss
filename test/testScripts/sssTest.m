classdef sssTest < matlab.unittest.TestCase
    
    properties 
        pwdPath
        sysCell
        deleteBenchmarks
        testPath
    end
    
    methods
        function testCase=sssTest()
            testCase.pwdPath=pwd;
            if exist('benchmarksSysCell.mat','file')
                testCase.deleteBenchmarks=0;
            else
                testCase.testPath=loadBenchmarks;
                testCase.deleteBenchmarks=1;
            end
            
            temp=load('benchmarksSysCell.mat');
            testCase.sysCell=temp.benchmarksSysCell;
            if isempty(testCase.sysCell)
                error('No benchmarks loaded.');
            end
            disp('sssTest executed...');
        end
    end
    
    methods(TestClassTeardown)
        function deleteBench(testCase)
            if testCase.deleteBenchmarks
                cd(testCase.testPath);
                delete('benchmarksSysCell.mat');
            end
            cd(testCase.pwdPath);
        end
    end
end