classdef testIterativeRefinement < sssTest
    
    methods (Test)  
        function testUsage(testCase)
            sys= sss('CDplayer');
            A = sys.A; B = sys.B; X = A\B;
            [L,U] = lu(A);
            
            % call default
            iterativeRefinement(A,B,X);
            [Xnew,rVec] = iterativeRefinement(A,B,X);
            
            verifySize(testCase, Xnew, size(X))
            verifyEqual(testCase, size(rVec,1), 1);
            
            % call with nonempty LU
            Opts.L = L; Opts.U = U;
            iterativeRefinement(A,B,X,Opts);
            
            % call with nonempty L but empty U
            clear Opts
            Opts.L = L;
            iterativeRefinement(A,B,X,Opts);
            
            % call with all sparse LU matrices
            clear Opts
            [L,U,P,Q,S] = lu(A);
            Opts.L = L; Opts.U = U; Opts.P = P; Opts.Q = Q; Opts.S = S;
            iterativeRefinement(A,B,X,Opts);

            
        end
        function testVectorBNoLu(testCase)
           % load model
           sys = sss('CDplayer');
           sys = sys(1,1); %make siso -> B is a vector
           
           % solve lse
           A = sys.A; B = sys.B;
           X = A\B; 
                      
           % introduce perturbation to degrade accuracy
           X = X + 1e6*eps*X;
           R = A*X-B; rNormPert = norm(R,'fro')/norm(B,'fro');
           
           % WILKINSON 
           Opts.method = 'wilkinson'; Opts.tol = 1e-15;
           [X,rNormVec] = iterativeRefinement(A,B,X,Opts);
           R = A*X-B; rNormRef = norm(R,'fro')/norm(B,'fro');
           
           verifyEqual(testCase, rNormRef, rNormVec(end),'RelTol',1e-10)
           verifyLessThanOrEqual(testCase, rNormRef, rNormPert);
           verifyLessThanOrEqual(testCase, rNormRef, Opts.tol); 
           
           % CGS 
           Opts.method = 'cgs'; Opts.tol = 1e-15;
           [X,rNormVec] = iterativeRefinement(A,B,X,Opts);
           R = A*X-B; rNormRef = norm(R,'fro')/norm(B,'fro');
           
           verifyEqual(testCase, rNormRef, rNormVec(end),'RelTol',1e-10)
           verifyLessThanOrEqual(testCase, rNormRef, rNormPert);
           verifyLessThanOrEqual(testCase, rNormRef, Opts.tol);
        end
        function testMatrixBNoLu(testCase)
           % load model
           sys = sss('CDplayer');
           
           % solve lse
           A = sys.A; B = sys.B;
           X = A\B; 
                      
           % introduce perturbation to degrade accuracy
           X = X + 1e6*eps*X;
           R = A*X-B; rNormPert = norm(R,'fro')/norm(B,'fro');
           
           % WILKINSON 
           Opts.method = 'wilkinson'; Opts.tol = 1e-15;
           [X,rNormVec] = iterativeRefinement(A,B,X,Opts);
           R = A*X-B; rNormRef = norm(R,'fro')/norm(B,'fro');
           
           verifyEqual(testCase, rNormRef, rNormVec(end),'RelTol',1e-10)
           verifyLessThanOrEqual(testCase, rNormRef, rNormPert);
           verifyLessThanOrEqual(testCase, rNormRef, Opts.tol); 
           
           % CGS 
           Opts.method = 'cgs'; Opts.tol = 1e-15;
           [X,rNormVec] = iterativeRefinement(A,B,X,Opts);
           R = A*X-B; rNormRef = norm(R,'fro')/norm(B,'fro');
           
           verifyEqual(testCase, rNormRef, rNormVec(end),'AbsTol',1e-10)
           verifyLessThanOrEqual(testCase, rNormRef, rNormPert);
           verifyLessThanOrEqual(testCase, rNormRef, Opts.tol);
        end
        function testMatrixBLu(testCase)
           % load model
           sys = sss('CDplayer');
           
           % solve lse
           A = sys.A; B = sys.B;
           [L,U,P,Q,S] = lu(A);
           X = A\B; 
                      
           % introduce perturbation to degrade accuracy
           X = X + 1e6*eps*X;
           R = A*X-B; rNormPert = norm(R,'fro')/norm(B,'fro');
           
           % WILKINSON 
           Opts = struct('method','wilkinson','tol',1e-15,...
                        'L',L, 'U', U, 'Q', Q, 'P', P, 'S', S);
           [X,rNormVec] = iterativeRefinement(A,B,X,Opts);
           R = A*X-B; rNormRef = norm(R,'fro')/norm(B,'fro');
           
           verifyEqual(testCase, rNormRef, rNormVec(end),'AbsTol',1e-10)
           verifyLessThanOrEqual(testCase, rNormRef, rNormPert);
           verifyLessThanOrEqual(testCase, rNormRef, Opts.tol); 
           
           % CGS 
           Opts.method = 'cgs';
           [X,rNormVec] = iterativeRefinement(A,B,X,Opts);
           R = A*X-B; rNormRef = norm(R,'fro')/norm(B,'fro');
           
           verifyEqual(testCase, rNormRef, rNormVec(end),'AbsTol',1e-10)
           verifyLessThanOrEqual(testCase, rNormRef, rNormPert);
           verifyLessThanOrEqual(testCase, rNormRef, Opts.tol);
        end
    end
end