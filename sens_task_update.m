% Study measurement upatre of sensor tasking
% - pseudo
% - sensor contraints ----
% - probability of detection
clc
clear all
close all
Px=[10,8;8,10];
mux=[5;5];
fn=2;


alphaa=pi;
Rmax=4;
dirn=0;
FOV=[alphaa,Rmax,dirn];
xsenspos=[4;5];
R=[5,-3;-3,5];
hn=2;
% [y,G,tp]= hxrb(xtarg,xsenspos,FOV);


%% Kalman filter update
[X,w]=UT_sigmapoints(mux,Px,2);

Y=zeros(size(X,1),hn);
G=zeros(size(X,1),1);
for i=1:1:size(X,1)
    [y,g,tp]= hxy(X(i,:)',xsenspos,FOV);
    Y(i,:)=y;
    G(i)=g;
end
RR=0;
for i=1:1:size(X,1)
    r=sqrt(sum((xsenspos-X(i,:)').^2));
    g= max((r/Rmax)^5,1)
    RR=RR+w(i)*g*R;
    %     if G(i)==-1
%         RR=RR+w(i);
%     else
%         RR=RR+w(i);
%     end
end
[R,RR]
% RR=R;
% [R,RR]

[mz,Pz]=MeanCov(Y,w);
Pz=Pz+RR;
Pcc=CrossCov(X,mux,Y,mz,w);
[xk,Pk]=KalmanUpdate(mux,Px,mz,Pz,Pcc,-1);
% [mk1,Pk1]=KalmanUpdate_modf_angle(mk,Pk,mz,Pz,Pcc,ym,ytp)


%%
figure(1)
plot_1sig_ellipse(mux,Px,'b')
hold on
plot_1sig_ellipse(mux,Pk,'r')
plot(X(:,1),X(:,2),'go','linewidth',2)
plot_circle_sensor(xsenspos(1),xsenspos(2),Rmax,'k')
% axis([-5,10,-5,10])
grid on

% N=100000;
% Xmc=mvnrnd(mux,Px,N);
% Y=zeros(N,2);
% G=zeros(N,1);
% for j=1:1:N
%     [y,g,tp]= hxy(Xmc(j,:)',xsenspos,FOV);
%     Y(j,:) =y;
%     G(j)=g;
% end
% Yin=Y(G==1,:);
% figure(2)
% hist3(Yin,[50,50])
% xlabel('r')
% ylabel('th')




