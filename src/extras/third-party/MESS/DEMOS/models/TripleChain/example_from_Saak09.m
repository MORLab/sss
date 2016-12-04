%%
n1=150;
alpha=0.01;
beta=alpha;
v=5e0;

%%
[M,D,K]=triplechain_MSD(n1,alpha,beta,v);
B=ones(3*n1+1,1);
Cp=B';
Cv=zeros(size(Cp));
%Cv=B';

%%
nsample=200;
%w=logspace(-4,0,nsample);
w=logspace(-4,2,nsample);

tro=zeros(1,nsample);
fprintf(['Computing TFMs of original systems and ' ...
  'MOR errors\n'])

%%
for k=1:nsample
  fprintf('\r Step %3d / %3d',k,nsample)
    Go = (Cp + 1i*w(k)*Cv)/(-w(k)*w(k)*M + 1i*w(k)*D + K) * B;
    tro(k) = max(svds(Go));
end
fprintf('\n\n');
figure(1)
loglog(w, tro)
xlabel('\omega')
ylabel('\sigma_{max}(G(j\omega))')
title('Transfer functions of original systems')
