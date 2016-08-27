function XY=sigmaEllipse(mu,P,n)
i=1;
XY=zeros(361,2);
for th=0:1:360
    XY(i,:)=(n*sqrtm(P)*[cosd(th);sind(th)]+mu(:))';
    i=i+1;
end

end