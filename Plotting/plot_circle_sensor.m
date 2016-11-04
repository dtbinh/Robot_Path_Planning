function plot_circle_sensor(xc,yc,R,c)

k=1;
X=zeros(length(0:0.1:360),2);
for th=0:0.1:360
    X(k,:)=[xc;yc]+R*[cosd(th);sind(th)];
    k=k+1;
end
% muv=A*[1;0]+mux;
plot(xc,yc,'k^',X(:,1),X(:,2),c,'linewidth',2,'MarkerSize',10)
% Xpt(:,1),Xpt(:,2),'go',
% Xpt [mux(1),muv(1)+mux(1)],[mux(2),muv(2)+mux(2)]
end