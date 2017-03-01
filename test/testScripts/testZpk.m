classdef testZpk < sssTest
    % testZpk - testing of zpk.m
   
    methods(Test)
        function testZpkObject(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};

                % call zpk for the sys
                kP = 8;
                kZ = 8;
                zpkData = zpk(sys,kP,'sm',kZ,'la');
                
                % check class
                if ~isa(zpkData,'zpk')
                   error('Wrong class!. Object should be of class "zpk"') 
                end
                
                % check correct number of zeros and poles
                for l = 1:sys.p
                    for j = 1:sys.m                        
                        % check correct number of zeros   
                        if length(zpkData(l,j).z{1,1}) ~= kZ
                            error('Wrong number of zeros!');
                        end
                        
                        % check correct number of poles                        
                        if length(zpkData(l,j).p{1,1}) ~= kP
                            error('Wrong number of poles!');
                        end
                    end
                end  
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-6,'AbsTol',1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end