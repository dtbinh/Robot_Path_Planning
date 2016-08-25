function plot_2Darc(xc,yc,alphaa,R,th,c,meth)
% th is the rotation matrix
% c is the colour 
% xc,yc is the vertex of cone
% equilateral triangle
% alpha is half angle, a is the center range 
% as the triag is plotted is upwar
% th=th-pi;
% R=[cos(th),-sin(th);sin(th),cos(th)];
% x=[0,a*tan(alphaa),-a*tan(alphaa),0];
% y=[0,a,a,0];
% % keyboard
% for i=1:1:length(x)
%      xy=R*[x(i);y(i)]+[xc;yc];
% x(i)=xy(1);
% y(i)=xy(2);
% end
% keyboard
theta0=th-alphaa;
theta1=th+alphaa;
theta=linspace(theta0,theta1,10);

x=[xc,zeros(1,length(theta))];
y=[yc,zeros(1,length(theta))];

for i=1:1:length(theta)
   x(i+1)= xc+R*cos(theta(i));
    y(i+1)= yc+R*sin(theta(i));
end
x=[x,xc];
y=[y,yc];

plot(x,y,c,'linewidth',2)
if strcmpi(meth,'fill')
    fill(x,y,c)
    alpha 0.2;
end
