function template_gen(x0,d,umin,umax,Vmin,Vmax,T,Vdx,Vdy)
close all
warning off
tic
%d is the grid length
%T is the time interval
%Vdx and Vdy are the disturbances in the x and y directions
% always assumed to start at the origin
% ending will be at the reachable grid nodes

% options=optimset('disp','off','TolCon',1e-7,'MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-5,'TolX',1e-5,'Jacobian','on');
options=optimset('TolCon',1e-7,'MaxFunEvals',1e5,'MaxIter',1e5,'TolFun',1e-5,'TolX',1e-5);


N=10; % variables to solve for i.e. N regualr spaces interval in [0,T]
dt=T/(N-1);

%XD=[[d,d,pi/2];[d;d;0]';[0;d;pi/2]';[0;d;0]';[d;0;pi/2]';[d;0;0]';...
%    [0;2*d;pi/2]';[0;2*d;0]';[2*d;0;pi/2]';[2*d;0;0]';[2*d;0;-pi/2]';[1*d;0;-pi/2]';...
%    [0;0;0]';[0;0;-pi/2]';[d;d;-pi/2]';[d;d;pi]';[d;0;pi]';[2*d;d;pi/2]'];

XD=[[d,d,pi/2];[d;d;0]';[0;d;pi/2]';[d;0;0]';[d;0;-pi/2]';[-d,d,pi];[-d,d,pi/2];[-d,0,pi];[-d,0,3*pi/2]];

%% initialize the control vector
templates1=cell(size(XD,1),1);
templates2=cell(size(XD,1),1);
templates3=cell(size(XD,1),1);
templates4=cell(size(XD,1),1);

%XD=vertcat(XD,repmat([-1,1,1],size(XD,1),1).*XD);

%%  xd=[d;d;pi/2]; 
parfor i=1:size(XD,1)
    i
    NN=10000;
U1=linspace(umin,umax,NN);
u1=U1(randi(NN,1,N));
U2=linspace(Vmin,Vmax,NN);
u2=U2(randi(NN,1,N));

u0=[u1,u2]';
lb=[umin*ones(N,1);Vmin*ones(N,1)];
ub=[umax*ones(N,1);Vmax*ones(N,1)];

    para=0;
    para.x0=x0;
para.d=d;
para.ulim=[umin,umax];
para.Vlim=[Vmin,Vmax];
para.T=T;
para.Vdist=[Vdx,Vdy];
para.J=0;
Rw=@(t)5*t*diag([5,1]);
Sw=diag([50,50,5000]);
xd=XD(i,:)';
[u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0,Vmin),u0,[],[],[],[],lb,ub,[]);


u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
[tJ,Jj]=ode45(@(t,x)[u1(t),(u2(t)-Vmin)]*Rw(t)*[u1(t);(u2(t)-Vmin)],[0,T],0);
J=Jj(end);
para.xd=xd;para.J=J;
templates1{i}=para;
templates2{i}=x;
templates3{i}=u;
templates4{i}=linspace(0,T,N);

end


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

% %%  xd=[d;d;0]; 
% Rw=[1:1:N]'/100;
% Sw=diag([10,10,1000]);
% xd=[d;d;0];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{2,1}=para;
% templates{2,2}=x;
% templates{2,3}=u;
% templates{2,4}=linspace(0,T,N);


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[0;d;pi/2]; 
% Rw=[1:1:N]'/100;
% Sw=diag([10,10,1000]);
% xd=[0;d;pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{3,1}=para;
% templates{3,2}=x;
% templates{3,3}=u;
% templates{3,4}=linspace(0,T,N);
% 


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[0;d;0]; 
% Rw=[1:1:N]'/100;
% Sw=diag([10,10,1000]);
% xd=[0;d;0];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{4,1}=para;
% templates{4,2}=x;
% templates{4,3}=u;
% templates{4,4}=linspace(0,T,N);
% 


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[d;0;pi/2]; 
% Rw=[1:1:N]'/100;
% Sw=diag([20,20,2000]);
% xd=[d;0;pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{5,1}=para;
% templates{5,2}=x;
% templates{5,3}=u;
% templates{5,4}=linspace(0,T,N);



%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[d;0;0]; 
% Rw=[1:1:N]'/100;
% Sw=diag([10,10,1000]);
% xd=[d;0;0];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{6,1}=para;
% templates{6,2}=x;
% templates{6,3}=u;
% templates{6,4}=linspace(0,T,N);
% 


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[0;2*d;pi/2]; 
% Rw=[1:1:N]'/100;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[0;2*d;pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{7,1}=para;
% templates{7,2}=x;
% templates{7,3}=u;
% templates{7,4}=linspace(0,T,N);


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[0;2*d;0]; 
% Rw=[1:1:N]'/100;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([20,20,2000]);
% xd=[0;2*d;0];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{8,1}=para;
% templates{8,2}=x;
% templates{8,3}=u;
% templates{8,4}=linspace(0,T,N);
% 
% 

%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[2*d;0;pi/2]; 
% Rw=[1:1:N]'/100;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[2*d;0;pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{9,1}=para;
% templates{9,2}=x;
% templates{9,3}=u;
% templates{9,4}=linspace(0,T,N);
% 
% 

%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[2*d;0,0]; 
% Rw=[1:1:N]'/100;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[2*d;0;0];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{10,1}=para;
% templates{10,2}=x;
% templates{10,3}=u;
% templates{10,4}=linspace(0,T,N);
% 
% 

%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[2*d;0;-pi/2]; 
% Rw=[1:1:N]'/10;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[2*d;0;-pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{11,1}=para;
% templates{11,2}=x;
% templates{11,3}=u;
% templates{11,4}=linspace(0,T,N);



%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[1*d;0;-pi/2]; 
% Rw=[1:1:N]'/10;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[d;0;-pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{12,1}=para;
% templates{12,2}=x;
% templates{12,3}=u;
% templates{12,4}=linspace(0,T,N);



%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[0;0;0]; 
% Rw=[1:1:N]'/1000;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[0;0;0];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{13,1}=para;
% templates{13,2}=x;
% templates{13,3}=u;
% templates{13,4}=linspace(0,T,N);
% 


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[0;0;-pi/2]; 
% Rw=[1:1:N]'/100;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[0;0;-pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{14,1}=para;
% templates{14,2}=x;
% templates{14,3}=u;
% templates{14,4}=linspace(0,T,N);
% 


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[d;d;-pi/2]; 
% Rw=[1:1:N]'/100;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[d;d;-pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{15,1}=para;
% templates{15,2}=x;
% templates{15,3}=u;
% templates{15,4}=linspace(0,T,N);



%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)


%%  xd=[d;d;pi]; 
% Rw=[1:1:N]'/1000;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[d;d;pi];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{16,1}=para;
% templates{16,2}=x;
% templates{16,3}=u;
% templates{16,4}=linspace(0,T,N);
% 


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[d;0;pi]; 
% Rw=[1:1:N]'/100;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[d;0;-pi];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{17,1}=para;
% templates{17,2}=x;
% templates{17,3}=u;
% templates{17,4}=linspace(0,T,N);
% 

%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)

%%  xd=[2*d;d;pi/2]; 
% Rw=[1:1:N]'/100;
% % Rw=@(t)diag([1,1])/10;
% Sw=diag([50,50,5000]);
% xd=[2*d;d;pi/2];
% [u,J]=fmincon(@(u)dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0),u0,[],[],[],[],lb,ub,[],options);
% 
% 
% u1 = @(t)interp1(linspace(0,T,N),u(1:N)',t,'cubic');
% u2 = @(t)interp1(linspace(0,T,N),u(N+1:end)',t,'cubic');
% [t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2,Vdx,Vdy),[0,T],x0);
% para.xd=xd;para.J=J;
% templates{18,1}=para;
% templates{18,2}=x;
% templates{18,3}=u;
% templates{18,4}=linspace(0,T,N);


%figure
%plot(x(:,1),x(:,2),'r--',x0(1),x0(2),'ks',xd(1),xd(2),'bo','linewidth',2,'MarkerSize',6)
toc
disp('done lol')
ddd=strcat(num2str(x0(1)),'_',num2str(x0(2)),'_',num2str(x0(3)),'s',num2str(d),'s',num2str(umin),'s',num2str(umax),'s',num2str(Vmin),'s',num2str(Vmax),'s',num2str(T),'s',num2str(Vdx),'s',num2str(Vdy));
ddd(findstr(ddd,'.'))='p';
save(ddd,'templates1','templates2','templates3','templates4','XD')
end
