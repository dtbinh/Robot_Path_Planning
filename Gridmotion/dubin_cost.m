function J=dubin_cost(u,T,Vdx,Vdy,Rw,xd,Sw,x0,Vmin)
% xd is the desired final destination as soft constraint
N=length(u)/2; %number of control variables..i.e time discretizations
u1=u(1:N);
u2=u(N+1:1:end);
U1 = @(t)interp1(linspace(0,T,N),u1,t,'cubic');
U2 = @(t)interp1(linspace(0,T,N),u2,t,'cubic');

% calculate control cost
% J1=sum((u1.^2).*Rw)+sum((u2.^2).*Rw);
%Rw=@(t)5*t*diag([5,1]);
[t,j]=ode45(@(t,x)[U1(t),(U2(t)-Vmin)]*Rw(t)*[U1(t);(U2(t)-Vmin)],[0,T],0);
J1=j(end);
% tt=linspace(0,T,N)';
% J1=sum((u1.^2).*Rw.*tt)+sum((u2.^2).*Rw.*tt);

% calculate the final time state cost

[t,x]=ode45(@(t,x)cont_dubin(t,x,U1,U2,Vdx,Vdy),[0,T],x0);
J2=(x(end,:)-xd')*Sw*(x(end,:)-xd')';

J=J1+J2;
end



