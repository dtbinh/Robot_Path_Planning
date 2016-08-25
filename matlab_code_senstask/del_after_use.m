function template_gen(x0,d,umin,umax,Vmin,Vmax,T,Vdx,Vdy)
%d is the grid length
%T is the time interval
%Vdx and Vdy are the disturbances in the x and y directions
% always assumed to start at the origin
% ending will be at the reachable grid nodes

% options=optimset('disp','iter','TolCon',1e-7,'MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-5,'TolX',1e-5,'Jacobian','on');
options=optimset('disp','iter','TolCon',1e-7,'MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-5,'TolX',1e-5);
%% initialize the control vector
templates=cell(1,4);
N=10; % variables to solve for i.e. N regualr spaces interval in [0,T]
dt=T/(N-1);

NN=10000;
U1=linspace(umin,umax,NN);
u1=U1(randi(NN,1,N));
U2=linspace(Vmin,Vmax,NN);
u2=U2(randi(NN,1,N));

u0=[u1,u2]';
lb=[umin*ones(N,1);Vmin*ones(N,1)];
ub=[umax*ones(N,1);Vmax*ones(N,1)];


%%  xd=[d;d;pi/2]; 
Rw=@(t)diag([1,1])/100;
Sw=diag([10,10,1000]);
xd=[d;d;pi/2];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{1,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{1,2}=x;
templates{1,3}=u;
templates{1,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[d;d;0]; 
Rw=@(t)diag([1,1])/100;
Sw=diag([10,10,1000]);
xd=[d;d;0];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);

u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{2,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{2,2}=x;
templates{2,3}=u;
templates{2,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[0;d;pi/2]; 
Rw=@(t)diag([1,1])/100;
Sw=diag([10,10,1000]);
xd=[0;d;pi/2];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{3,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{3,2}=x;
templates{3,3}=u;
templates{3,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[0;d;0]; 
Rw=@(t)diag([1,1])/100;
Sw=diag([10,10,1000]);
xd=[0;d;0];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{4,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{4,2}=x;
templates{4,3}=u;
templates{4,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[d;0;pi/2]; 
Rw=@(t)diag([1,1])/100;
Sw=diag([20,20,2000]);
xd=[d;0;pi/2];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{5,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{5,2}=x;
templates{5,3}=u;
templates{5,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[d;0;0]; 
Rw=@(t)diag([1,1])/100;
Sw=diag([10,10,1000]);
xd=[d;0;0];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{6,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{6,2}=x;
templates{6,3}=u;
templates{6,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[0;2*d;pi/2]; 
Rw=[1:1:N]'/100;
% Rw=@(t)diag([1,1])/10;
Sw=diag([50,50,5000]);
xd=[0;2*d;pi/2];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{7,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{7,2}=x;
templates{7,3}=u;
templates{7,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[0;2*d;0]; 
Rw=[1:1:N]'/100;
% Rw=@(t)diag([1,1])/10;
Sw=diag([20,20,2000]);
xd=[0;2*d;0];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{8,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{8,2}=x;
templates{8,3}=u;
templates{8,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[2*d;0;pi/2]; 
Rw=[1:1:N]'/100;
% Rw=@(t)diag([1,1])/10;
Sw=diag([50,50,5000]);
xd=[2*d;0;pi/2];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{9,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{9,2}=x;
templates{9,3}=u;
templates{9,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[2*d;0,0]; 
Rw=[1:1:N]'/100;
% Rw=@(t)diag([1,1])/10;
Sw=diag([50,50,5000]);
xd=[2*d;0;0];
u=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
templates{10,1}=[x0',d,umin,umax,Vmin,Vmax,T,Vdx,Vdy];
templates{10,2}=x;
templates{10,3}=u;
templates{10,4}=linspace(0,T,N);

figure
plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


