function [y,x_,m,k,index,L,U,p] = simInit(A,B,C,D,E,u,x,Ts,Ts_sample,isDescriptor,oneOutput)

if oneOutput
    x_      = [];
    m       = inf;
    k       = []; 
    index   = [];
else
    m       = round(Ts_sample/Ts);
    x_      = zeros(length(A),round(size(u,1)/m));
    x_(:,1) = x; %initial state
    k       = 1; %index for state
    index   = 1;
end

% Precompute LU decomposition of E
if isDescriptor
    [L,U,p] = lu(E,'vector');
else
    L = []; U = []; p = [];
end

% Initialize output
y       = zeros(size(C,1),size(u,1));
y(:,1)  = C*x + D*u(1,:)';