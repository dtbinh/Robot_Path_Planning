% template_gen([0;0;pi/2],30,-1,1,10,15,5,0,0)
% template_gen([0;0;pi/2],50,-1,1,10,15,5,0,0)
close all
% load('0_0_1p5708s30s-1s1s10s15s5s0s0')
II=1:size(templates1,1);
for i=1:size(templates1,1)
%     templates{18,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
JJ=templates1{II(i),1}.J;
x0=templates1{II(i),1}.x0;
xd=templates1{II(i),1}.xd;
x=templates2{II(i),1};
u=templates3{II(i),1};
TT=templates4{II(i),1};
tt=linspace(0,TT(end),200);
% templates{18,3}=u;
% templates{18,4}=linspace(0,T,N);

% figure(i)

% subplot(2,1,1)
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo',[x0(1),x0(1)],[x0(2),x0(2)+5],'k',[xd(1),xd(1)+5*cos(xd(3))],[xd(2),xd(2)+5*sin(xd(3))],'b','linewidth',3,'MarkerSize',8)
hold on
plot(XD(:,1),XD(:,2),'bo',-XD(:,1),XD(:,2),'bo','MarkerSize',5,'linewidth',2)
axis([-50,50,-10,50])
set(gca,'XTick',-40:40:40);
set(gca,'YTick',[-40,0,40]);
axis square
plot_prop_paper
grid on
% pause(0.5)
% saveas(gcf,strcat(num2str(i),'d40'),'fig')
% saveas(gcf,strcat(num2str(i),'d40'),'pdf')
% pause(0.5)
% subplot(2,1,2)
% u1 = @(t)interp1(TT,u(1:length(u)/2)',t,'cubic');
% u2 = @(t)interp1(TT,u(length(u)/2+1:end)',t,'cubic');
% plot(tt,u1(tt),tt,u2(tt))
% legend('u','V')
% title(num2str(JJ))
% if i>5
% axis([-20,20,-5,45])
% else
% axis([-5,45,-5,45])
% end

% xlabel('x')
% ylabel('y')
% plot_prop_paper

% figure(i+6)
% % subplot(2,6,i+6)
% plot(-x(:,1),x(:,2),'r--',-x0(1),x0(2),'ks',-xd(1),xd(2),'bo',[x0(1),x0(1)],[x0(2),x0(2)+5],'k',[-xd(1),-xd(1)-5*cos(-xd(3))],[xd(2),xd(2)-5*sin(-xd(3))],'b','linewidth',2,'MarkerSize',6)
% hold on
% plot(XD(:,1),XD(:,2),'bo',-XD(:,1),XD(:,2),'bo','MarkerSize',5,'linewidth',2)
% if i==4 || i==3
% axis([-20,20,-5,45])
% else
% axis([-45,5,-5,45])
% end
% axis square
% grid on
% plot_prop_paper
end
% 

%% generate the lookup table


