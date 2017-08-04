function [y,x_,k,index] = simUpdate(y,x_,k,index,x,u,i,m,C,D)    

% Update output and state
y(:,i) = C*x + D*u(i,:)';
if ~isempty(x_)
    if mod(i,m) == 0
        k       = k+1;            
        x_(:,k) = x;
        index   = [index i];
    end
end