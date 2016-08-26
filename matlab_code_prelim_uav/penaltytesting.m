[phi,alpha]=meshgrid(-90*pi/180:1*pi/180:(90)*pi/180,2*pi/180:0.5*pi/180:10*pi/180);
Lr=0;
Lt=0;
sens.xc1=[0,0];
for i=1:1:size(phi,1)
    for j=1:1:size(phi,2)
LL=sensor2penalty([2,2],[phi(i,j),alpha(i,j)],sens);
Lr(i,j)=LL(1,1);
Lt(i,j)=LL(2,2);
    end
end
mesh(phi*180/pi,alpha*180/pi,Lt)
xlabel('phi')
ylabel('alpha')