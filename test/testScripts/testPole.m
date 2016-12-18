classdef testPole < sssTest
    % testPoles - testing of pole.m
   
    methods(Test)
        function testLM(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};

                if ~sys.isDae
                    Opts.type='lm';
                    k=20;
                    p=pole(sys,k,Opts);     
                    
                    % compare results of eig and eigs
                    pEig=eig(full(sys.A),full(sys.E));

                    % remove poles at infinity
                    pEig=pEig(abs(real(pEig))<1e6);

                    tbl=table(-abs(p),p);
                    tbl=sortrows(tbl);
                    p=tbl.p;
                    tbl=table(-abs(pEig),pEig);
                    tbl=sortrows(tbl);
                    pEig=tbl.pEig;

                    %remove single complex element (for comparison)
                    if abs(imag(sum(p)))>1e-12
                       p(abs(imag(p)-imag(sum(p)))<1e-16)=[];
                    end
                    
                    % make sure sorting is correct
                    p=cplxpair(p);
                    pEig=cplxpair(pEig(1:size(p,1)));
                    
                    verification(testCase, p,pEig(1:size(p,1)));
                end
            end
        end
        function testSM(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if ~sys.isDae
                    Opts.type='sm';
                    k=10;
                    p=pole(sys,k,Opts);

                    % compare results of eig and eigs
                    pEig=eig(full(sys.A),full(sys.E));

                    % remove poles at infinity
                    pEig=pEig(abs(real(pEig))<1e6);

                    tbl=table(abs(p),p);
                    tbl=sortrows(tbl);
                    p=tbl.p;
                    tbl=table(abs(pEig),pEig);
                    tbl=sortrows(tbl);
                    pEig=tbl.pEig;

                    % remove single complex element (for comparison)
                    if abs(imag(sum(p)))>1e-12
                       p(abs(imag(p)-imag(sum(p)))<1e-16)=[];
                    end

                    verification(testCase, p,pEig(1:size(p,1)));
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-6,'AbsTol',1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end