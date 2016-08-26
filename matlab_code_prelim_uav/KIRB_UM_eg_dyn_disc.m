function xk1=KIRB_UM_eg_dyn_disc(xk,para)
T=para(1);
% T=1;
% keyboard
xk1=[1,0,T,0;
     0,1,0,T;
     0,0,1,0;
     0,0,0,1]*xk;
end