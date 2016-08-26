%% testing the gaussian sum mixture evaluation of expected KL divergence
% two radial radars
clear
P1=[6,5;5,6];
 P2=[6,-5;-5,10];
 mux=[35,30;30,35];
% mux=[5,3];
Px=[reshape(P1,1,2^2);reshape(P2,1,2^2)];
wx=[0.5;1-0.5];
sigma = cat(3,P1,P2);
 obj = gmdistribution(mux,sigma,wx');
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
xc2=[10,0];

model.R=R;
model.fn=2;
model.hn=2;
model.hx=@(x,para)[sqrt((x(1)-xc1(1))^2+(x(2)-xc1(2))^2);sqrt((x(1)-xc2(1))^2+(x(2)-xc2(2))^2)];
model.gx=@(x,para)[sensor_radial_penalty(x,xc1,para(1:2)),0;0,sensor_radial_penalty(x,xc2,para(3:4))];

ang=0:5*pi/180:1*pi/2;

i=1;
GMM.mu=mux;
GMM.P=Px;
GMM.w=wx;
DKL_gsm=0;
DKL_gsmod=0;
alp=[2*pi/180,5*pi/180,10*pi/180];
DKL_gsmodp=0;
for phi1=ang
      for alpha1=alp
        for phi2=ang
            for alpha2=alp

DKL_gsmod=ExptKL_modf_MeasOPT_gsm(GMM,model,[phi1,alpha1,phi2,alpha2]);
if DKL_gsmod>DKL_gsmodp
    para=[phi1,alpha1,phi2,alpha2];
DKL_gsmodp=DKL_gsmod;
end
i=i+1
            end
        end
      end
      
 end

[x,y]=meshgrid(20:1:40);
CC=0;
for i=1:1:length(x)
    for j=1:1:length(x)
        CC(i,j)=pdf(obj,[x(i,j) y(i,j)]);
    end
end
contour(x,y,CC,25)
plot_sens_view(xc1,para(2),para(1),80,'b')
plot_sens_view(xc2,para(4),para(3),80,'r')
axis([0,50,0,50])

