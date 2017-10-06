classdef testSecondToFirst < sssTest
    
    methods (Test)  
        function testDifferentOptions(testCase)
            % test the differnt options for the F-matrix
            load('LF10');
            
            % option 'I'
            Opts.transf2nd = 'I';
            sys_I = second2first(M,D,K,B,C,[],Opts);
            
            % option 'K'
            Opts.transf2nd = 'K';
            sys_K = second2first(M,D,K,B,C,[],Opts);
            
            % option '-K'
            Opts.transf2nd = '-K';
            sys_minK = second2first(M,D,K,B,C,[],Opts);
            
            % optoin 'alpha'
            Opts.transf2nd = 'alpha';
            sys_alpha = second2first(M,D,K,B,C,[],Opts);
            
            % option scalar value
            Opts.transf2nd = 5;
            sys_scalar = second2first(M,D,K,B,C,[],Opts);
            
            % option user defined matrix
            F = rand(18);
            F = (F+F');
            F = F + 20*eye(18);
            Opts.transf2nd = F;
            sys_matrix = second2first(M,D,K,B,C,[],Opts);
            
            % create bode-plots
            [mag_I,phase_I,omega] = bode(sys_I);
            [mag_K,phase_K] = bode(sys_K,omega);
            [mag_minK,phase_minK] = bode(sys_minK,omega);
            [mag_alpha,phase_alpha] = bode(sys_alpha,omega);
            [mag_scalar,phase_scalar] = bode(sys_scalar,omega);
            [mag_matrix,phase_matrix] = bode(sys_matrix,omega);
            
            % check if the bode-plots are equal
            verifyEqual(testCase,mag_I,mag_K,'RelTol',1e-8, ...
                'Verification failed for option "K" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_K,'RelTol',1e-8, ...
                'Verification failed for option "K" (phase)!');
            verifyEqual(testCase,mag_I,mag_minK,'RelTol',1e-8, ...
                'Verification failed for option "-K" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_minK,'RelTol',1e-8, ...
                'Verification failed for option "-K" (phase)!');
            verifyEqual(testCase,mag_I,mag_alpha,'RelTol',1e-8, ...
                'Verification failed for option "alpha" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_alpha,'RelTol',1e-8, ...
                'Verification failed for option "alpha" (phase)!');
            verifyEqual(testCase,mag_I,mag_scalar,'RelTol',1e-8, ...
                'Verification failed for option "scalar" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_scalar,'RelTol',1e-8, ...
                'Verification failed for option "scalar" (phase)!');
            verifyEqual(testCase,mag_I,mag_matrix,'RelTol',1e-8, ...
                'Verification failed for option "matrix" (magnitude)!');
            verifyEqual(testCase,phase_I,phase_matrix,'RelTol',1e-8, ...
                'Verification failed for option "matrix" (phase)!');     
        end
        
        function testMatrixOutput(testCase)
            load('gyro')
            sys         =   second2first(M,1e-6*K,K,B,C);
            [A,B,C,D,E] =   second2first(M,1e-6*K,K,B,C);
            sys2 = sss(A,B,C,D,E);
            verifyClass(testCase, sys,'sss');
            verifyEqual(testCase, sys, sys2);
        end
    end
end