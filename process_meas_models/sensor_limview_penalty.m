function L=sensor_limview_penalty(x,xc,Rmax,para)
% for now alpha is fixed
phi=para(1);% look direction
alpha=para(2);
% 
% Rmid=Rmax/2;
 r=sqrt((x(1)-xc(1))^2+(x(2)-xc(2))^2);
% d=abs(Rmid-r);
% % 
if r>Rmax
    L=100;
    return;
end

beta=atan2(x(2)-xc(2),x(1)-xc(1));
y=abs(phi-beta);

if y>alpha
    L=100;
    return; 
end

L=1+(1*y/alpha)^1;



end