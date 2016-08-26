%% testing the gaussian sum mixture evaluation of expected KL divergence
% single sensor with range and bearing with 2 cost functions.
clear
clc
% P1=[10,5;5,10];
  P2=[1,0;0,15];
%  mux=[35,30;20,35];
 mux=[35,35];
 Px=reshape(P2,1,2^2);
 wx=1;
sigma = cat(3,P2);
 obj = gmdistribution(mux,sigma,wx');
 
% GMM.mu=[35,30];
% GMM.P=reshape(P2,1,4);
% GMM.w=1;
GMM.mu=mux;
for i=1:1:length(wx)
GMM.P(i,:)=reshape(obj.Sigma(:,:,i),1,4);
end
GMM.w=obj.PComponents';
GMM=GMM_moms(GMM);

%  [x,y]=meshgrid(-10:0.5:10);
%  for i=1:1:length(x)
%      for j=1:1:length(y)
%          Y(i,j)=pdf(obj,[x(i,j) y(i,j)]);
%      end
%  end
%  surf(x,y,Y)
%  
%  ezsurf(@(x,y)pdf(obj,[x y]),[0 40],[0 40])

R=[2,0;0,2];

xc1=[0,0];
model.lam=0;

model.R=1;
model.fn=2;
model.hn=1;
% model.hx=@(x,para)[sqrt((x(1)-xc1(1))^2+(x(2)-xc1(2))^2);atan2(x(2)-xc1(2),x(1)-xc1(1))];
% model.gx=@(x,para)[sensor_radial_penalty(x,xc1,para),0;0,1];
model.hx=@(x,para)atan2(x(2)-xc1(2),x(1)-xc1(1));
model.gx=@(x,para)sensor_radial_penalty(x,xc1,para);
model.marg_nx=0;
ang=0:2*pi/180:1*pi/2;

i=1;

%GMM2.mu=GMM.mean';
%GMM2.P=reshape(GMM.cov,1,4);
%GMM2.w=1;

DKL_gsm=0;
DKL_gsmod=0;
alp=[2*pi/180,5*pi/180,10*pi/180,20*pi/180];
DKL_gsmodp=0;
DKL_1gcp=0;
D=0;
i=0;
j=0;
for phi1=ang
i=i+1;
j=0;
     for alpha1=alp
j=j+1;

DKL_gsmod=ExptKL_modf_MeasOPT_gsm(GMM,model,[phi1,alpha1]);
D(i,j)=DKL_gsmod;
% DKL_1gc=ExptKL_modf_MeasOPT_gsm(GMM2,model,[phi1,alpha1]);
if DKL_gsmod>DKL_gsmodp
    para=[phi1,alpha1];
DKL_gsmodp=DKL_gsmod;
end

% if DKL_1gc>DKL_1gcp
%     para_1gcp=[phi1,alpha1];
% DKL_1gcp=DKL_1gc;
% end



      end
      
end
%  figure
% plot(ang,D(:,1))
% figure
[x,y]=meshgrid(10:1:50);
CC=0;
CC2=0;
for i=1:1:length(x)
    for j=1:1:length(x)
        CC(i,j)=pdf(obj,[x(i,j) y(i,j)]);
%         CC2(i,j)=mvnpdf([x(i,j),y(i,j)],GMM.mean',GMM.cov);
    end
end
contour(x,y,CC,25)
% hold on
% contour(x,y,CC2,5,'--')
% hold off
plot_sens_view(xc1,para(2),para(1),80,'b')
% plot_sens_view(xc1,para_1gcp(2),para_1gcp(1),80,'y')
axis([0,50,0,50])

[DKL_gsmodp,para*180/pi]