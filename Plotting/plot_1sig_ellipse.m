function plot_1sig_ellipse(mux,Px,c)
%  [Xpt,w]=GH_points(mux,Px,2);

A=sqrtm(Px);
k=1;
for th=0:0.1:360
    X=[cosd(th);sind(th)];
    Y(k,:)=A*X+mux;
    k=k+1;
end
% muv=A*[1;0]+mux;
plot(Y(:,1),Y(:,2),c,'linewidth',2)
% Xpt(:,1),Xpt(:,2),'go',
% Xpt [mux(1),muv(1)+mux(1)],[mux(2),muv(2)+mux(2)]
end