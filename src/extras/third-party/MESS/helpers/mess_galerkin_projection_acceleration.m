function Z=mess_galerkin_projection_acceleration(Z,type,eqn,oper,fopts)
%  Galerkin projection onto subspace spanned by low-rank factor Z of ADI
%  method for solving FXE^T + EXF^T = -GG^T.
%
% Input:
%  Z                 Low-rank factor Z
%
%  eqn               structure with data for A, E, G
%                          eqn.E(optional, eqn.haveE specifies whether it is
%                          there) in the above equation with ZZ' approximating X
%
%  oper              structure contains function handles for
%                          operations with A, E
%
%  opts              options structure that should contain following members
%  opts.ortho        implicit or explitic orthogonalization of Z defaults to 1,
%                          i.e., explicit orthogonalization via orth().
%  opts.meth         method for solving projected Lyapunov equation
%
%
%  opts.comprType    compression type, default is 0 for SVD based
%                          compression, to use RRQR based compression
%
%  fopts
%
% Output:
%  Z                 Updated solution factor Z after prolongation
%
% uses oparatorfunctions mul_A, mul_E

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) Jens Saak, Martin Koehler, Peter Benner and others 
%               2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
%

factorize=1;

switch type
    case {'lyapunov','lyap','lyapc','lyapchol'}
        opts=fopts.adi;
    case {'riccati','care','care_nwt_fac'}
        opts=fopts.nm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts.projection,'ortho') || isempty(opts.projection.ortho)
    opts.projection.ortho=1;
end
if ~isfield(opts.projection,'meth') || isempty(opts.projection.meth)
    switch type
        case {'lyapunov','lyap','lyap2solve'}
            opts.projection.meth='lyap2solve';
        case {'lyapc','lyapchol'}
            opts.projection.meth='lyapchol';
        case {'riccati','care'}
            opts.projection.meth='care';
        case {'care_nwt_fac' }
            opts.projection.meth='care_nwt_fac';
        case {'mess_dense_nm'}
            opts.projection.meth='mess_dense_nm';
    end
end

if ~isfield(fopts,'rosenbrock'), fopts.rosenbrock=[]; end
if isstruct(fopts.rosenbrock)&&isfield(fopts.rosenbrock,'tau')
    rosenbrock = 1;
    pc = -1 / (2 * fopts.rosenbrock.tau);
else
    rosenbrock = 0;
end
if ~isfield(fopts,'bdf'), fopts.bdf=[]; end
if isstruct(fopts.bdf) && isfield(fopts.bdf, 'tau') && isfield(fopts.bdf, 'beta')
    bdf = 1;
    pc = -1 / (2 * fopts.bdf.tau * fopts.bdf.beta);
else
    bdf = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute projector matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.projection.ortho
    Z=orth(Z);
else
    [U,S,~]=svd(full(Z'*Z));
    s=diag(S);
    sk=find(s>eps*s(1), 1, 'last' );
    Z=Z*U(:,1:sk)*diag(1./sqrt(s(1:sk)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the projected Matrix equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch type
    case {'lyapunov','lyap','lyapchol','lyap_sgn_fac'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The Lyapunov equation case
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lyapunov = 1;
        % Project down matrices
        if eqn.type=='N'
            if bdf
                A=Z'*((fopts.bdf.tau * fopts.bdf.beta) * ...
                    oper.mul_ApE(eqn, fopts,'N',pc, 'N', Z,'N'));
            elseif rosenbrock
                A=Z'*(oper.mul_ApE(eqn, fopts,'N',pc, 'N', Z,'N') ...
                    - eqn.U * (eqn.V' * Z));
            else
                A=Z'*(oper.mul_A(eqn, fopts,'N',Z,'N'));
            end
            B=Z'*eqn.B;
        elseif eqn.haveUV
            if bdf
                A=Z'*((fopts.bdf.tau * fopts.bdf.beta) * ...
                    oper.mul_ApE(eqn, fopts,'T',pc, 'T',Z,'N')-eqn.V*(eqn.U'*Z));
            elseif rosenbrock
                A=Z'*(oper.mul_ApE(eqn, fopts,'T',pc, 'T',Z,'N')-eqn.V*(eqn.U'*Z));
            else
            A=Z'*(oper.mul_A(eqn, fopts,'T',Z,'N')-eqn.V*(eqn.U'*Z));
            end
            B=Z'*eqn.G;
        else
            if bdf
                A=Z'*((fopts.bdf.tau * fopts.bdf.beta) * ...
                    oper.mul_ApE(eqn, fopts,'T',pc, 'T',Z,'N'));
            elseif rosenbrock
                A=Z'*(oper.mul_ApE(eqn, fopts,'T',pc, 'T',Z,'N') ...
                    - eqn.V * (eqn.U' * Z));
            else
                A=Z'*(oper.mul_A(eqn, fopts,'T',Z,'N'));
            end
            B=Z'*eqn.C';
        end
        
        if eqn.haveE
            if eqn.type=='N'
                E=Z'*(oper.mul_E(eqn, fopts,'N',Z,'N'));
            else
                E=Z'*(oper.mul_E(eqn, fopts,'T',Z,'N'));
            end
        else
            E=[];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Choose solver for the small equation
        file_not_found = 0;
        switch opts.projection.meth
            case 'lyapchol'
                if exist('lyapchol','file')
                    XC=lyapchol(A,B,E);
                    factorize=0;
                else
                    file_not_found=1;
                end
                
            case 'lyap_sgn_fac'
                if exist('lyap_sgn_fac','file')
                    XC=lyap_sgn_fac(A,B',E);
                    factorize=0;
                else
                    file_not_found=1;
                end
                
            case {'lyap','lyapunov'}
                if exist('lyap','file')
                    X=lyap(A,B*B',[],E);
                else
                    file_not_found=1;
                end
                
            case 'lyap2solve'
                if exist('lyap2solve','file')
                    X=lyap2solve(A,B*B',[],E);
                else
                    file_not_found=1;
                end
            otherwise
                file_not_found=1;
        end
        if file_not_found
            warning('MESS:GP_missing_solver',...
                ['Galerkin projection acceleration skipped due to missing direct'...
                ' Lyapunov solver, or unknown method']);
            factorize=0;
        end
    case {'riccati','care','care_nwt_fac','mess_dense_nm'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The Riccati equation case
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        file_not_found = 0;
        lyapunov = 0;
        
        if bdf
            A=Z'*((fopts.bdf.tau * fopts.bdf.beta) * ...
                oper.mul_ApE(eqn, fopts,'N',pc, 'N',Z,'N'));
        elseif rosenbrock
            A=Z'*(oper.mul_ApE(eqn, fopts,eqn.type,pc, 'N',Z,'N') ...
                - eqn.U * (eqn.V' * Z));
        else
            if eqn.setUV
                A = Z' * (oper.mul_A( eqn, fopts, 'N', Z, 'N' ) ...
                    - eqn.U(:,1:eqn.UVsize) * (eqn.V(:,1:eqn.UVsize)' * Z));
            else
                A=Z'*(oper.mul_A(eqn, fopts,'N',Z,'N'));
            end
        end
        B=Z'*eqn.B;
        C=eqn.C*Z;
        if eqn.haveE
            M=Z'*(oper.mul_E(eqn, fopts,'N',Z,'N'));
        else
            M=[];
        end
        
        switch opts.projection.meth
            case {'care', 'riccati'}
                if exist('care','file')
                    if ~isempty(M)
                        X=care(A,B,C'*C,eye(size(B,2)),[],M);
                    else
                        X=care(A,B,C'*C,eye(size(B,2)));
                    end
                else
                    file_not_found=1;
                end
                
            case 'care_nwt_fac'
                if exist('care_nwt_fac','file')
                    if ~isempty(M)
                        XC = care_nwt_fac([],M\A,M\B,C,1e-12,50);
                        XC = XC*M;
                    else
                        XC = care_nwt_fac([],A,B,C,1e-12,50);
                    end
                else
                    file_not_found=1;
                end
                factorize=0;
                
            case 'mess_dense_nm'
                if exist('mess_dense_nm','file')
                    if ~isempty(M)
                        X=mess_dense_nm(A,B,C,M);
                    else
                        X=mess_dense_nm(A,B,C,[]);
                    end
                else
                    file_not_found=1;
                end
            otherwise
                file_not_found=1;
        end
        if file_not_found
            warning('MESS:GP_missing_solver',...
                ['Galerkin projection acceleration skipped due to missing'...
                ' Riccati solver']);
            factorize=0;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the projected solution was not already computed in factored form
% compute a symmetric factorization now and update the large factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if factorize
    if (exist('cholp','file'))
        [XC,P,I]=cholp(X);
        XC=XC*P';
        if I && lyapunov
            warning('MESS:proj_sol_semidef',...
                'The solution of the projected equation was semidefinite.')
        end
    else
        [~,S,V]=svd(X);
        s=diag(S);
        r=find(s>s(1)*eps);
        XC=diag(sqrt(s(r)))*V(:,r)';
    end
end
if exist('XC','var'), Z=Z*XC'; end
