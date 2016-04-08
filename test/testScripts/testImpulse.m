classdef testImpulse < sssTest
    
    methods (Test)
        function testImpulse1(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                t=0:0.1:5;
                
                [actH,actT]=impulse(sys,t);
                [expH,expT]=impulse(ss(sys),t);
                
                
                actSolution={actH,actT};
                expSolution={expH,expT};
                
                verification (testCase, actSolution, expSolution);
                verifyInstanceOf(testCase, actH , 'double', 'Instances not matching');
                verifyInstanceOf(testCase, actT , 'double', 'Instances not matching');
                verifySize(testCase, actH, size(expH), 'Size not matching');
                verifySize(testCase, actT, size(expT), 'Size not matching');
                
            end
        end
%         function testImpulsePlot(testCase)
%             load('CDplayer.mat')
%             sys=sss(A,B,C);
%             t=0:0.01:5;
%             sys.InputName = {'u1','u2'}';
%             sys.OutputName = {'y1','y2'}';
%             sysSS = ss(sys);
%             sys1_1 = sys(1,1);
%             sys1_2 = sys(1,2);
%             sys1_12 = sys(1,:);
%             impulse(sys,sysSS,sys1_1,sys1_2,sys1_12,t);
%             title('testImpulsePlotE')
%         end
%         function testImpulsePlotE(testCase)
%             load('CDplayer.mat')
%             sys=sss(2*A,2*B,C,[],2*eye(size(A)));
%             sys.InputName = {'u1','u2'}';
%             sys.OutputName = {'y1','y2'}';
%             sysSS = ss(sys);
%             sys1_1 = sys(1,1);
%             sys1_2 = sys(1,2);
%             sys1_12 = sys(1,:);
%             figure
%             impulse(sys,sysSS,sys1_1,sys1_2,sys1_12);
%             title('testImpulsePlotE')
%         end
    end
end

function [] = verification (testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
    'Difference between actual and expected exceeds relative tolerance');
end