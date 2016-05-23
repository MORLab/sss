classdef testZeros < sssTest
    % testZeros - testing of zeros.m
   
    methods(Test)
        function testLM(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if sys.isSiso && ~sys.isDae

    %                 sys=loadSss('heat-cont');
                    disp(sys.Name);
                    Opts.type='lm';
                    k=20;


                    z=zeros(sys,k,Opts);

                    % compare results of eig and eigs
                    zEig=eig(full([sys.A,sys.B;sys.C,sys.D]),[full(sys.E),zeros(sys.n,sys.m);zeros(sys.p,sys.n),zeros(sys.p,sys.m)]);

                    % remove zeros at infinity
                    zEig=zEig(abs(real(zEig))<1e6);

                    if strcmp(Opts.type,'lm')
                        tbl=table(-abs(z),z);
                        tbl=sortrows(tbl);
                        z=tbl.z;
                        tbl=table(-abs(zEig),zEig);
                        tbl=sortrows(tbl);
                        zEig=tbl.zEig;
                    else
                        tbl=table(abs(z),z);
                        tbl=sortrows(tbl);
                        z=tbl.z;
                        tbl=table(abs(zEig),zEig);
                        tbl=sortrows(tbl);
                        zEig=tbl.zEig;
                    end

                    % remove single complex element (for comparison)
                    if abs(imag(sum(z)))>1e-12
                       z(abs(imag(z)-imag(sum(z)))<1e-16)=[];
                    end

                    verification(testCase, z,zEig(1:size(z,1)));
                end
            end
        end
        function testSM(testCase)
            for i=1:length(testCase.sysCell)
                sys=testCase.sysCell{i};
                if sys.isSiso && ~sys.isDae

                    disp(sys.Name);
                    Opts.type='sm';
                    k=20;


                    z=zeros(sys,k,Opts);

                    % compare results of eig and eigs
                    zEig=eig(full([sys.A,sys.B;sys.C,sys.D]),[full(sys.E),zeros(sys.n,sys.m);zeros(sys.p,sys.n),zeros(sys.p,sys.m)]);

                    % remove zeros at infinity
                    zEig=zEig(abs(real(zEig))<1e6);

                    if strcmp(Opts.type,'lm')
                        tbl=table(-abs(z),z);
                        tbl=sortrows(tbl);
                        z=tbl.z;
                        tbl=table(-abs(zEig),zEig);
                        tbl=sortrows(tbl);
                        zEig=tbl.zEig;
                    else
                        tbl=table(abs(z),z);
                        tbl=sortrows(tbl);
                        z=tbl.z;
                        tbl=table(abs(zEig),zEig);
                        tbl=sortrows(tbl);
                        zEig=tbl.zEig;
                    end

                    % remove single complex element (for comparison)
                    if abs(imag(sum(z)))>1e-12
                       z(abs(imag(z)-imag(sum(z)))<1e-16)=[];
                    end
                    
                    verification(testCase, z,zEig(1:size(z,1)));
                end
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution,  expSolution,'RelTol',1e-6,'AbsTol',1e-6,...
    'Difference between actual and expected exceeds relative tolerance');
end