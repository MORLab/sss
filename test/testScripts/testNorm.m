classdef testNorm < sssTest
    % testNorm - testing of norm.m
    %
    % Description:
    %   The function norm.m is tested (3 tests) on:
    %    + Norm of a SISO benchmark system.
    %    + Norm of a SISO random system.
    %    + Norm of a MISO random system.
    %    + Norm of a SIMO random system.
    %    + Norm of MIMO benchmark system.
    %    + Norm of a MIMO random system.
    %    + Verifies for every case the following inputs/outputs (norm(sys),
    %       norm(sys,inf), [n,fpeak]=norm(sys,inf), norm(sys,2)
    %
    % ------------------------------------------------------------------
    %   This file is part of sssMOR, a Sparse State Space, Model Order
    %   Reduction and System Analysis Toolbox developed at the Institute
    %   of Automatic Control, Technische Universitaet Muenchen.
    %   For updates and further information please visit www.rt.mw.tum.de
    %   For any suggestions, submission and/or bug reports, mail us at
    %                     -> sssMOR@rt.mw.tum.de <-
    % ------------------------------------------------------------------
    % Authors:      Alessandro Castagnotto
    %               Jorge Luiz Moreira Silva
    % Last Change:  26 Out 2015
    % Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
    % ------------------------------------------------------------------
    
    methods(Test)
        function testNorm1(testCase)
            for i=1:length(testCase.sysCell)
                sys_sss=testCase.sysCell{i};
                sys_ss=ss(sys_sss);
                actNorm1=norm(sys_sss);
                [actNorm2,actFreq]=norm(sys_sss,inf);
                sys_sss.D=zeros(size(sys_sss.D));
                actNorm3=norm(sys_sss,2);
                actNorm = [actNorm1,actNorm2,actFreq,actNorm3];
                expNorm1=norm(sys_ss);
                [expNorm2,expFreq]=norm(sys_ss,inf);
                sys_ss.D=zeros(size(sys_ss.D));
                expNorm3=norm(sys_ss,2);
                if norm(freqresp(sys_ss,actFreq),2)>expNorm2
                    expNorm2=actNorm2;
                    expFreq=actFreq;
                end
                if (actFreq==0 &&expFreq==0)
                    expNorm2=actNorm2;
                end
                expNorm = [expNorm1,expNorm2,expFreq,expNorm3];
                verification(testCase, actNorm, expNorm);
            end
        end
    end
end

function [] = verification(testCase, actSolution, expSolution)
verifyEqual(testCase, actSolution(1:4),  expSolution(1:4),'RelTol',1e-3,...
    'Difference between actual and expected exceeds relative tolerance');
end