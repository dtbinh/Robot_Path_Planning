clear 
% clc
i=1; 
P1=[2,1;1,20];
m1=[5;3];
ang=0:0.3*pi/180:2*pi;
% [X,W]=GH_points(m1,P1,7);
% [X,W]=UT_sigmapoints(m1,P1,2);
[X,W]=conjugate_dir_gausspts(m1,P1);
for phi=ang
D1(i)=expt_KL_sens(phi,m1,P1,X,W);    
D2(i)=expt_KL_sens2(phi,m1,P1,X,W);
i=i+1;
end
plot(ang*180/pi,D1,ang*180/pi,D2)
legend('1','2')
[a1,b1]=max(D1);
[a2,b2]=max(D2);
[ang(b1),ang(b2),atan2(m1(2),m1(1))]*180/pi