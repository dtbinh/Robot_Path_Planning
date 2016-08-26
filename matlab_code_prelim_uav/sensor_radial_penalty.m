function L=sensor_radial_penalty(x,xc,para)
phi=para(1);% look direction
alpha=para(2); % wide angle look

beta=atan2(x(2)-xc(2),x(1)-xc(1));
y=abs(phi-beta);
if y<=alpha
    g=50+(alpha)*(y)*500;
else
    g=1e3;
    
end
L=g;
end