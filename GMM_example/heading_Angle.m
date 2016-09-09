clear 
clc
time.t0 = 0;
time.dt = 5;
time.tf = 125+90+125+30+125;
time.tspan = time.t0 : time.dt : time.tf;
time.nSteps = length(time.tspan);
%__________________________________________________________________________

%% ------------------------------------------------------------------------
% model

T=time.dt;
model.para_dt(1)=T;
model.fx = @(x)CT_turn_model(x,model.para_dt);
model.x0tr=[25000,10000,-120,0,0,pi]';
x0=[25000,10000,-120,0,0,pi]';
i=1;
x_mc=zeros(100,6);
x_mc(i,:)=x0;
for tt=0:T:495
    i=i+1;
     if tt<125
         x0(5)=0;
     elseif tt>=125 && tt<125+90
         x0(5)=1*pi/180;
     elseif tt>=125+90 && tt<125+90+125
         x0(5)=0;
     elseif tt>=125+90+125 && tt<125+90+125+30
         x0(5)=-3*pi/180;
     elseif tt>=125+90+125+30 && tt<125+90+125+125
         x0(5)=0;
     end
    [t,x0]=ode45_discc(@model.fx,tt,T,tt+5,x0,1e-200);
x_mc(i,:)=x0(end,:);
x0=x0(end,:)';
    tt
    plot(x_mc(1:i,1),x_mc(1:i,2))
axis([-1.7,2.8,-1.7,1.5]*10^4)

end
