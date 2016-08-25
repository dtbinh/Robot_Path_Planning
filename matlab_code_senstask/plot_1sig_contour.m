function plot_1sig_contour(mux,Px)
%  [Xpt,w]=GH_points(mux,Px,2);
 [Xpt,w]=UT_sigmapoints(mux,Px,2);
Y=zeros(length(0:0.1:360),2);
muv=mux(3:4);
muv=muv/norm(muv);

mux=mux(1:2);
mux=mux(:);
Px=Px(1:2,1:2);

A=sqrtm(Px);
k=1;
for th=0:0.1:360
    X=[cosd(th);sind(th)];
    Y(k,:)=A*X+mux;
    k=k+1;
end
% muv=A*[1;0]+mux;
plot(mux(1),mux(2),'ks',Y(:,1),Y(:,2),'g','linewidth',2)
% Xpt(:,1),Xpt(:,2),'go',
% Xpt [mux(1),muv(1)+mux(1)],[mux(2),muv(2)+mux(2)]
end