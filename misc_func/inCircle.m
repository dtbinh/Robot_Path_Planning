function Y=inCircle(X,xc,yc,R)
% check if the points X are within the cricle, 
% return  1 for yes and 0 for no
% The points in X
n=size(X,1);
Y=zeros(n,1);
% keyboard
S=(X(:,1)-xc).^2+(X(:,2)-yc).^2-R^2;
Y(S>0)=0;
Y(S<=0)=1;

end
