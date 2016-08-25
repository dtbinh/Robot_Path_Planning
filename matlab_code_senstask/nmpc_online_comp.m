u1=@(t)sin(0.3*t);
u2=@(t)[5];
[t,x]=ode45(@(t,x)cont_dubin(t,x,u1,u2),[0:0.02:20],[0,0,0]);

figure(1)
plot(x(:,1),x(:,2))
figure(2)
plot(t,u1(t))

%% Solve the complete NMPC iteratively here



%% Linearize the model and solve the NMPC