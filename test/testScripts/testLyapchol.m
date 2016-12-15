classdef testLyapchol < sssTest
    methods(Test)
        function testLyapchol1(testCase) 
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    [S,R]=lyapchol(sys);
                    verification(testCase, sys, S, R);
                end
            end
        end
        function testLyapcholOpts(testCase)
            Opts.method='adi';
            Opts.q=120;
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae && sys.n>300
                    [S,R]=lyapchol(sys,Opts);
                    verification(testCase, sys, S, R);
                end
            end
        end
        function testLyapcholE(testCase) 
            for i=1:2
                if i<2
                    sys=loadSss('LF10'); % E is symmetric
                else
                    sys=loadSss('SpiralInductorPeec');
                end
                [S,R]=lyapchol(sys);               
                verification(testCase, sys, S, R);      
            end
        end
    end
end

function [] = verification(testCase,sys,  S, R)
    X = S*S'; Y = R*R';
    tol = 1e-3;
    verifyLessThan(testCase, ...
        norm(sys.A*(X)*sys.E' + sys.E*(X)*sys.A' + sys.B*sys.B'),...
        tol);
    verifyLessThan(testCase, ...
        norm(sys.A'*(Y)*sys.E + sys.E'*(Y)*sys.A + sys.C'*sys.C),...
        tol);

end