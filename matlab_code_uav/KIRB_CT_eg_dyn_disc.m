function xk1=KIRB_CT_eg_dyn_disc(xk,T)
% global T
% T=para(1);
% T=1;
xk=xk(:);
omg=xk(end,1)+1e-10;
% keyboard
xk1=[1,0,sin(omg*T)/omg,-(1-cos(omg*T))/omg,0;...
     0,1,(1-cos(omg*T))/omg,sin(omg*T)/omg,0;
     0,0,cos(omg*T),-sin(omg*T),0;
     0,0,sin(omg*T),cos(omg*T),0;
     0,0,0,0,1]*xk;
% xk1=xk1';
end
%% discrete dynamic equations for MC
