function L=sensor_radial_penalty2(x,xc,Rmax,para)
% phi=para(1);% look direction
% alpha=para(2); % wide angle look

Rmid=Rmax/2;
r=sqrt((x(1)-xc(1))^2+(x(2)-xc(2))^2);
 d=abs(Rmid-r);

if r>Rmax
    L=1e4;
    return;
end

L=d^2;
end