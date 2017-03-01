classdef testFreqresp < sssTest
    
    methods (Test)  
        function testFreqresp1(testCase) %real frequencies
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                if ~sys_sss.isDae
                    sys_ss=ss(sys_sss);
                    t=0:100:1000;
                    t = [-t,t];

                    [actG, actOmega]=freqresp(sys_sss,t);
                    [expG, expOmega]=freqresp(sys_ss,t);

                    actSolution={actG, actOmega};
                    expSolution={expG, expOmega};

                    verification (testCase, actSolution, expSolution);
                    verifyInstanceOf(testCase, actG , 'double', 'Instances not matching');
                    verifyInstanceOf(testCase, actOmega , 'double', 'Instances not matching');
                    verifySize(testCase, actG, size(expG), 'Size not matching');
                    verifySize(testCase, actOmega, size(expOmega), 'Size not matching');
                end
            end
        end
        function testFreqresp2(testCase) %imaginary frequencies
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                if ~sys_sss.isDae
                    sys_ss=ss(sys_sss);
                    t=1:100:1000;
                    t=[-1i*t,1i*t];

                    [actG, actOmega]=freqresp(sys_sss,t);
                    [expG, expOmega]=freqresp(ss(sys_ss),t);


                    actSolution={actG, actOmega};
                    expSolution={expG, expOmega};

                    verification (testCase, actSolution, expSolution);
                    verifyInstanceOf(testCase, actG , 'double', 'Instances not matching');
                    verifyInstanceOf(testCase, actOmega , 'double', 'Instances not matching');
                    verifySize(testCase, actG, size(expG), 'Size not matching');
                    verifySize(testCase, actOmega, size(expOmega), 'Size not matching');
                end
            end
        end
        function testFreqresp3(testCase) %complex frequencies
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                if ~sys_sss.isDae
                    sys_ss=ss(sys_sss);
                    t=1e2*randn(1,10) + 1e4*1i*randn(1,10);

                    [actG, actOmega]=freqresp(sys_sss,t);
                    [expG, expOmega]=freqresp(ss(sys_ss),t);


                    actSolution={actG, actOmega};
                    expSolution={expG, expOmega};

                    verification (testCase, actSolution, expSolution);
                    verifyInstanceOf(testCase, actG , 'double', 'Instances not matching');
                    verifyInstanceOf(testCase, actOmega , 'double', 'Instances not matching');
                    verifySize(testCase, actG, size(expG), 'Size not matching');
                    verifySize(testCase, actOmega, size(expOmega), 'Size not matching');
                end
            end
        end
        function testOmegaRange(testCase)
            % input of frequnecy range
            
            for i=1:length(testCase.sysCell)
                sys  = testCase.sysCell{i};
                if ~sys.isDae
                    w = {10,100};

                    [~,omega] = freqresp(sys,w);
                    verifyEqual(testCase,omega(1),w{1},  'Wrong frequency returned');
                    verifyEqual(testCase,omega(end),w{2},'Wrong frequency returned');
                end
            end
        end
        function testStaticGain(testCase)
            % test a static gain model
            sys1=sss([],[],[],5);
            [actSolution1,omega]=freqresp(sys1);
            expSolution1=freqresp(ss(sys1),omega);
            
            sys2=sss(-speye(3),zeros(3,2),zeros(3),[1,2;3,4;5,6]);
            [actSolution2,omega]=freqresp(sys2);
            expSolution2=freqresp(ss(sys2),omega);
            verification(testCase,{actSolution1,actSolution2},{expSolution1,expSolution2});
        end
    end
end
    
function [] = verification (testCase, actSolution, expSolution)
          verifyEqual(testCase, actSolution, expSolution, 'RelTol', 0.1,'AbsTol',0.000001, ...
               'Difference between actual and expected exceeds relative tolerance');
end