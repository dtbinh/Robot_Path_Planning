function R=sensor_fov_R(x,para,xc,Rmax)

phi=para(1);% look direction
alpha=para(2);


beta=atan2(x(2)-xc(2),x(1)-xc(1));
y=abs(phi-beta);

if y>alpha
    La=1000;

else

La=1+(y/alpha)^1;
end

Rmid=Rmax/2;
r=sqrt((x(1)-xc(1))^2+(x(2)-xc(2))^2);
 d=abs(Rmid-r);

if r>Rmax
    Lr=1000;

else

Lr=1+(d/Rmid)^1;
end

if La>=900 || Lr>=900
 R=[100,0;0,100];
else
R=[Lr,0;0,La];
end