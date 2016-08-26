function f=Lpen(xc,x,phi,alpha)
B=atan2(x(2)-xc(2),x(1)-xc(1));
y=abs(B-phi);
% if y<=alpha
%     f=2*(y)^2;
% else
%     f=10*(y-alpha)+2*alpha^2;
% end
f=10*(y)^2;
end

