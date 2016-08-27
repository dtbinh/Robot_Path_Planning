function plot_ellipse_cov(mu,P,n)
% n is the sigma variance multiples at which ellipses are drawn
A=sqrtm(P);
X=zeros(length(n),length(0:0.1:360));
Y=zeros(length(n),length(0:0.1:360));
k=1;
for th=0:0.1:360
    
    B=repmat(mu,length(n),1)+repmat(n(:),1,length(mu)).*repmat((A*[cosd(th);sind(th)])',length(n),1);
    X(:,k)=B(:,1);
    Y(:,k)=B(:,2);
    k=k+1;
end
plot(X',Y','linewidth',2)