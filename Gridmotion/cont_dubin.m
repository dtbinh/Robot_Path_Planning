function dx=cont_dubin(t,x,u1,u2,Vdx,Vdy)
% u(1) turn rate control
% u(2) velocity control
% if x(3)>pi
%     x(3)=x(3)-pi;
% end
% if x(3)<-pi
%     x(3)=x(3)+pi;
% end

dx=zeros(3,1);
dx(1)=u2(t)*cos(x(3))+Vdx;
dx(2)=u2(t)*sin(x(3))+Vdy;
dx(3)=u1(t);
end
