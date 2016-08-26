%% function maximize divergence
%pdf at k+1
function D=expt_KL_sens(phi,m1,P1,X,W)
% P1=[4,2;2,6];
% m1=[10;10];
% Pxk1=@(x)1/((2*pi)^1*sqrt(abs(det(Pk1))))*exp(-0.5*(x-mk1)'*inv(Pk1)*(x-mk1));

%measuremetn model
sigr=1;
sigth=1*pi/180;
R=[sigr^2,0;0,sigth^2];
L=@(phi,x)(1+(phi-atan2(x(2),x(1)))^2);
% R=@(x,phi)[sigr^2,0;0,sigth^2]*(1+(phi-atan2(x(2),x(1)))^2);
h=@(x)[sqrt(x(1)^2+x(2)^2);atan2(x(2),x(1))];
% Pz=@(z,x,phi)1/((2*pi)^1*sqrt(abs(det(Rk1(x,phi)))))*exp(-0.5*(z-h(x))'*inv(Rk1(x,phi))*(z-h(x)));
%bayes rule at time k+1

Pz=0;
PL=0;
mz=0;
Pxz=0;
for i=1:1:length(X)
    mz=mz+W(i)*h(X(i,:));
end
for i=1:1:length(X)
    Pz=Pz+W(i)*(h(X(i,:))-mz)*(h(X(i,:))-mz)';
    PL=PL+W(i)*L(phi,X(i,:))^2;
    Pxz=Pxz+W(i)*(X(i,:)'-m1)*(h(X(i,:))-mz)';
end
Pz=Pz+PL*R;


K=Pxz*inv(Pz);
P2=P1-K*Pxz';
D=0.5*log(det(P2)/det(P1))-0.5*length(mz)+0.5*trace(inv(P2)*P1)+0.5*trace(K'*inv(P2)*K*Pz);
end