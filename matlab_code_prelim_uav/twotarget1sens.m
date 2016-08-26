%% testing the gaussian sum mixture evaluation of expected KL divergence
% single sensor with range and bearing -- 2 targets
clear
P=[0.001,0,0,0;0,0.001,0,0;0,0,10,0;0,0,0,10];

 GMM.mu=[35,25,25,35];
 GMM.P=reshape(P,1,16);
 GMM.w=1;
 % mux=[5,3];
% Px=[reshape(P1,1,2^2);reshape(P2,1,2^2)];
% wx=[0.6;1-0.6];
% sigma = cat(3,P1,P2);
%  obj = gmdistribution(mux,sigma,wx');
%  
% GMM.mu=obj.mu;
% for i=1:1:length(wx)
% GMM.P(i,:)=reshape(obj.Sigma(:,:,i),1,4);
% end
% GMM.w=obj.PComponents';
%  GMM=GMM_moms(GMM);
%  [x,y]=meshgrid(-10:0.5:10);
%  for i=1:1:length(x)
%      for j=1:1:length(y)
%          Y(i,j)=pdf(obj,[x(i,j) y(i,j)]);
%      end
%  end
%  surf(x,y,Y)
%  
%  ezsurf(@(x,y)pdf(obj,[x y]),[0 40],[0 40])

R=[1,0;0,1];

xc1=[0,0];
model.lam=0;

model.R=R;
model.fn=4;
model.hn=2;
model.hx=@(x,para)[atan2(x(2)-xc1(2),x(1)-xc1(1));atan2(x(4)-xc1(2),x(3)-xc1(1))];
model.gx=@(x,para)[sensor2penalty(x(1:2),xc1,para),0;0,sensor2penalty(x(3:4),xc1,para)];

ang=0:2*pi/180:1*pi/2;

i=1;

% GMM2.mu=GMM.mean';
% GMM2.P=reshape(GMM.cov,1,4);
% GMM2.w=1;

DKL_gsm=0;
DKL_gsmod=0;
alp=[2*pi/180,5*pi/180,10*pi/180];
DKL_gsmodp=0;
% DKL_1gcp=0;
for phi1=ang
      for alpha1=alp


DKL_gsmod=ExptKL_modf_MeasOPT_gsm(GMM,model,[phi1,alpha1]);
% DKL_1gc=ExptKL_modf_MeasOPT_gsm(GMM2,model,[phi1,alpha1]);
if DKL_gsmod>DKL_gsmodp
    para=[phi1,alpha1];
DKL_gsmodp=DKL_gsmod;
end
i=i+1
% if DKL_1gc>DKL_1gcp
%     para_1gcp=[phi1,alpha1];
% DKL_1gcp=DKL_1gc;
% end



      end
      
 end

[x,y]=meshgrid(10:1:50);
CC=0;
CC2=0;
P=reshape(GMM.P,4,4);
for i=1:1:length(x)
    for j=1:1:length(x)
%         CC(i,j)=pdf(obj,[x(i,j) y(i,j)]);
  CC(i,j)=mvnpdf([x(i,j),y(i,j)],GMM.mu(1:2),P(1:2,1:2));
        CC2(i,j)=mvnpdf([x(i,j),y(i,j)],GMM.mu(3:4),P(3:4,3:4));
    end
end
contour(x,y,CC,5)
hold on
contour(x,y,CC2,5,'--')
hold off
plot_sens_view(xc1,para(2),para(1),80,'b')
% plot_sens_view(xc1,para_1gcp(2),para_1gcp(1),80,'y')
axis([0,50,0,50])

para*180/pi