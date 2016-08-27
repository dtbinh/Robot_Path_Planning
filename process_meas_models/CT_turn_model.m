function xk1=CT_turn_model(xk,para,T)
% T=para(1);
% T=1;
omg=xk(5,1);
xk1=zeros(5,1);
if omg==0
    xk1(1:4,1)=[1,0,T,0;
     0,1,0,T;
     0,0,1,0;
     0,0,0,1]*xk(1:4,1);
 xk1(5,1)=omg;
%  xk1(6,1)=xk(6,1);
 return
end

xk1(1:5,1)=[1,0,sin(omg*T)/omg,-(1-cos(omg*T))/omg,0;...
     0,1,(1-cos(omg*T))/omg,sin(omg*T)/omg,0;
     0,0,cos(omg*T),-sin(omg*T),0;
     0,0,sin(omg*T),cos(omg*T),0;
     0,0,0,0,1]*xk(1:5,1);
%  xk1(6)=xk(6)+omg*T;
 
end