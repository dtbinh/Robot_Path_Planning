function L=sensor2penalty(x,xc,para)
phi=para(1);% look direction
alpha=para(2); % wide angle look

beta=atan2(x(2)-xc(2),x(1)-xc(1));
y=abs(phi-beta);
if y<=alpha
    g=2*alpha*y;
else
    g=(pi/4);
    
end
L=g;
end