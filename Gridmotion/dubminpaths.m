function L=dubminpaths(r0,rf,X0,Xf,p1,p2,dir)
L=2;
x0=X0(1);
y0=X0(2);
ph0=X0(3);
xf=Xf(1);
yf=Xf(2);
phf=Xf(3);

R=@(th)[cos(th),-sin(th);sin(th),cos(th)];
r1=r0;
r2=rf;


% lc at start
Xlc=[x0,y0]+(R(ph0+pi/2)*[r1;0])';
x1lc=Xlc(1);
y1lc=Xlc(2);

% rc circle at start 
Xrc=[x0,y0]+(R(ph0-pi/2)*[r1;0])';
x1rc=Xrc(1);
y1rc=Xrc(2);

% lc circlew at end
Xlc=[xf,yf]+(R(phf+pi/2)*[r2;0])';
x2lc=Xlc(1);
y2lc=Xlc(2);

% rc circle at end
Xrc=[xf,yf]+(R(phf-pi/2)*[r2;0])';
x2rc=Xrc(1);
y2rc=Xrc(2);


% d=sqrt((x0-xf)^2+(y0-yf)^2);
tgnt1=@(r1,r2,d)(sqrt(r1^2*d^2/(r1+r2)^2-r1^2)+sqrt(r2^2*d^2/(r1+r2)^2-r2^2));
tgnt2=@(r1,r2,d)sqrt(d^2-(r1-r2)^2);

%% dist rC1 to lC2
p1=sign(p1);
p2=sign(p2);
if p1>0 && p2<0 
xxc1=x1rc;yyc1=y1rc;
xxc2=x2lc;yyc2=y2lc;
end   
if p1>0 && p2>0 
xxc1=x1rc;yyc1=y1rc;
xxc2=x2rc;yyc2=y2rc;
end  
if p1<0 && p2<0 
xxc1=x1lc;yyc1=y1lc;
xxc2=x2lc;yyc2=y2lc;
end  
if p1<0 && p2>0 
xxc1=x1lc;yyc1=y1lc;
xxc2=x2rc;yyc2=y2rc;
end  
if p1>0 && p2<0 
xxc1=x1rc;yyc1=y1rc;
xxc2=x2lc;yyc2=y2lc;
end  

x1=0;
y1=0;
x2=xxc2-xxc1;
y2=yyc2-yyc1;
xt1=x0-xxc1;
yt1=y0-yyc1;
xt2=xf-xxc1;
yt2=yf-yyc1;
Alp1=atan2(yt1,xt1);
Alp2=atan2(yt2-y2,xt2-x2);
th=atan2(y2,x2);
Xx=R(-th)*[x2;y2];
x2=Xx(1);
y2=Xx(2);
ph1=ph0-th;
ph2=phf-th;
d1=norm([x2,y2]);
Alp1=Alp1-th;
Alp2=Alp2-th;

if dir==1
B1=atan(sqrt(r1^2*d1^2/(r1+r2)^2-r1^2)/r1);
if B1<0
    B1=B1+2*pi;
end
% line between the contacr points
x=(r1*cos(B1)):0.005:(d1+r2*cos(B1+pi));
y=(r2*sin(B1+pi))+(r1*sin(B1)-r2*sin(B1+pi))*(x-d1-r2*cos(B1+pi))/(r1*cos(B1)-d1-r2*cos(B1+pi));
XL=zeros(length(x),2);
for i=1:1:length(x)
   XL(i,:)=([xxc1;yyc1]+R(th)*[x(i);y(i)])'; 
end
inc_v=[x(end)-x(1);y(end)-y(2);0]/norm([x(end)-x(1);y(end)-y(2);0]);
r_v=[x(1)-x1;y(1)-y1;0]/norm([x(1)-x1;y(1)-y1;0]);
rot_dir1=-sum(sign(cross(r_v,inc_v)));
T2=norm([x(end)-x(1);y(end)-y(2);0]);
%plot 1st circles
% Alp1=ph1-pi/2;
if Alp1<0
    Alp1=Alp1+2*pi;
end
if Alp1>=B1 && rot_dir1==1 && p1==1
    tt=linspace(B1,Alp1,400);
elseif Alp1<B1 && rot_dir1==1 && p1==1
    tt=linspace(Alp1,Alp1-(2*pi-(B1-Alp1)),400);
elseif Alp1>=B1 && rot_dir1==-1 && p1==-1
    tt=linspace(B1,B1-(2*pi-(Alp1-B1)),400);
elseif Alp1<B1 && rot_dir1==-1 && p1==-1
    tt=linspace(Alp1,B1,400);    
else
     L=1e15;
    disp('not feasible')
    return
end
T1=r1*(max(tt)-min(tt));

XC1=zeros(400,2);
for i=1:1:400
   XC1(i,:)=[xxc1;yyc1]+R(th)*[r1*cos(tt(i));r1*sin(tt(i))];
end
% plot 2nd circles
% Alp=ph2-pi/2;
if Alp2<0
    Alp2=Alp2+2*pi;
end
B2=B1+pi;
if B2<0
    B2=B2+2*pi;
end
if B2>2*pi
    B2=B2-2*pi;
end

inc_v=[x(end)-x(1);y(end)-y(2);0]/norm([x(end)-x(1);y(end)-y(2);0]);
r_v=[x(end)-x2;y(end)-y2;0]/norm([x(end)-x2;y(end)-y2;0]);
rot_dir=-sum(sign(cross(r_v,inc_v)));

if Alp2>=B2 && rot_dir==1&& p2==1
    tt=linspace(B2,Alp2-2*pi,400);
elseif Alp2<B2 && rot_dir==1 && p2==1
    tt=linspace(B2,Alp2,400);
elseif Alp2>=B2 && rot_dir==-1 && p2==-1
     tt=linspace(B2,Alp2,400);
elseif Alp2<B2 &&  rot_dir==-1 && p2==-1
    tt=linspace(B2,2*pi+Alp2,400);         
else
    L=1e15;
    disp('not feasible')
    return
end
T3=r2*(max(tt)-min(tt));
XC2=zeros(400,2);
for i=1:1:400
   XC2(i,:)=[xxc2;yyc2]+R(th)*[r2*cos(tt(i));r2*sin(tt(i))];
end
L=T1+T2+T3;
end

if dir==2
B1=-atan(sqrt(r1^2*d1^2/(r1+r2)^2-r1^2)/r1);
if B1<0
    B1=B1+2*pi;
end
% line between the contacr points
x=(r1*cos(B1)):0.005:(d1+r2*cos(B1+pi));
y=(r2*sin(B1+pi))+(r1*sin(B1)-r2*sin(B1+pi))*(x-d1-r2*cos(B1+pi))/(r1*cos(B1)-d1-r2*cos(B1+pi));
XL=zeros(length(x),2);
for i=1:1:length(x)
   XL(i,:)=([xxc1;yyc1]+R(th)*[x(i);y(i)])'; 
end
T2=norm([x(end)-x(1);y(end)-y(2);0]);
inc_v=[x(end)-x(1);y(end)-y(2);0]/norm([x(end)-x(1);y(end)-y(2);0]);
r_v=[x(1)-x1;y(1)-y1;0]/norm([x(1)-x1;y(1)-y1;0]);
rot_dir1=-sum(sign(cross(r_v,inc_v)));

%plot 1st circles
% Alp1=ph1-pi/2;
if Alp1<0
    Alp1=Alp1+2*pi;
end
if Alp1>=B1 && rot_dir1==1 && p1==1
    tt=linspace(B1,Alp1,400);
elseif Alp1<B1 && rot_dir1==1 && p1==1
    tt=linspace(Alp1,Alp1-(2*pi-(B1-Alp1)),400);
elseif Alp1>=B1 && rot_dir1==-1 && p1==-1
    tt=linspace(B1,B1-(2*pi-(Alp1-B1)),400);
elseif Alp1<B1 && rot_dir1==-1 && p1==-1
    tt=linspace(Alp1,B1,400);    
else
     L=1e15;
    disp('not feasible')
    return
end
T1=r1*(max(tt)-min(tt));

XC1=zeros(400,2);
for i=1:1:400
   XC1(i,:)=[xxc1;yyc1]+R(th)*[r1*cos(tt(i));r1*sin(tt(i))];
end
% plot 2nd circles
% Alp=ph2-pi/2;
if Alp2<0
    Alp2=Alp2+2*pi;
end
B2=B1+pi;
if B2<0
    B2=B2+2*pi;
end
if B2>2*pi
    B2=B2-2*pi;
end
inc_v=[x(end)-x(1);y(end)-y(2);0]/norm([x(end)-x(1);y(end)-y(2);0]);
r_v=[x(end)-x2;y(end)-y2;0]/norm([x(end)-x2;y(end)-y2;0]);
rot_dir=-sum(sign(cross(r_v,inc_v)));

if Alp2>=B2 && rot_dir==1&& p2==1
    tt=linspace(B2,Alp2-2*pi,400);
elseif Alp2<B2 && rot_dir==1 && p2==1
    tt=linspace(B2,Alp2,400);
elseif Alp2>=B2 && rot_dir==-1 && p2==-1
     tt=linspace(B2,Alp2,400);
elseif Alp2<B2 &&  rot_dir==-1 && p2==-1
    tt=linspace(B2,2*pi+Alp2,400);         
else
    L=1e15;
    disp('not feasible')
    return
end
T3=r2*(max(tt)-min(tt));
XC2=zeros(400,2);
for i=1:1:400
   XC2(i,:)=[xxc2;yyc2]+R(th)*[r2*cos(tt(i));r2*sin(tt(i))];
end
L=T1+T2+T3;
end

if dir==3
B1=atan(sqrt(d1^2-(r1-r2)^2)/abs(r1-r2));
if B1<0
    B1=B1+2*pi;
end
% line between the contacr points
x=(r1*cos(B1)):0.005:(d1+r2*cos(B1));
y=(r2*sin(B1))+(r1*sin(B1)-r2*sin(B1))*(x-d1-r2*cos(B1))/(r1*cos(B1)-d1-r2*cos(B1));
XL=zeros(length(x),2);
for i=1:1:length(x)
   XL(i,:)=([xxc1;yyc1]+R(th)*[x(i);y(i)])'; 
end
T2=norm([x(end)-x(1);y(end)-y(2);0]);
inc_v=[x(end)-x(1);y(end)-y(2);0]/norm([x(end)-x(1);y(end)-y(2);0]);
r_v=[x(1)-x1;y(1)-y1;0]/norm([x(1)-x1;y(1)-y1;0]);
rot_dir1=-sum(sign(cross(r_v,inc_v)));

%plot 1st circles
% Alp1=ph1-pi/2;
if Alp1<0
    Alp1=Alp1+2*pi;
end
if Alp1>=B1 && rot_dir1==1 && p1==1
    tt=linspace(B1,Alp1,400);
elseif Alp1<B1 && rot_dir1==1 && p1==1
    tt=linspace(Alp1,Alp1-(2*pi-(B1-Alp1)),400);
elseif Alp1>=B1 && rot_dir1==-1 && p1==-1
    tt=linspace(B1,B1-(2*pi-(Alp1-B1)),400);
elseif Alp1<B1 && rot_dir1==-1 && p1==-1
    tt=linspace(Alp1,B1,400);    
else
     L=1e15;
    disp('not feasible')
    return
end
T1=r1*(max(tt)-min(tt));
XC1=zeros(400,2);
for i=1:1:400
   XC1(i,:)=[xxc1;yyc1]+R(th)*[r1*cos(tt(i));r1*sin(tt(i))];
end
% plot 2nd circles
% Alp=ph2-pi/2;
if Alp2<0
    Alp2=Alp2+2*pi;
end
B2=B1;
if B2<0
    B2=B2+2*pi;
end
if B2>2*pi
    B2=B2-2*pi;
end

inc_v=[x(end)-x(1);y(end)-y(2);0]/norm([x(end)-x(1);y(end)-y(2);0]);
r_v=[x(end)-x2;y(end)-y2;0]/norm([x(end)-x2;y(end)-y2;0]);
rot_dir=-sum(sign(cross(r_v,inc_v)));

if Alp2>=B2 && rot_dir==1&& p2==1
    tt=linspace(B2,Alp2-2*pi,400);
elseif Alp2<B2 && rot_dir==1 && p2==1
    tt=linspace(B2,Alp2,400);
elseif Alp2>=B2 && rot_dir==-1 && p2==-1
     tt=linspace(B2,Alp2,400);
elseif Alp2<B2 &&  rot_dir==-1 && p2==-1
    tt=linspace(B2,2*pi+Alp2,400);         
else
    L=1e15;
    disp('not feasible')
    return
end
T3=r2*(max(tt)-min(tt));
XC2=zeros(400,2);
for i=1:1:400
   XC2(i,:)=[xxc2;yyc2]+R(th)*[r2*cos(tt(i));r2*sin(tt(i))];
end
L=T1+T2+T3;
end

if dir==4
B1=-atan(sqrt(d1^2-(r1-r2)^2)/abs(r1-r2));
if B1<0
    B1=B1+2*pi;
end
% line between the contacr points
x=(r1*cos(B1)):0.005:(d1+r2*cos(B1));
y=(r2*sin(B1))+(r1*sin(B1)-r2*sin(B1))*(x-d1-r2*cos(B1))/(r1*cos(B1)-d1-r2*cos(B1));
XL=zeros(length(x),2);
for i=1:1:length(x)
   XL(i,:)=([xxc1;yyc1]+R(th)*[x(i);y(i)])'; 
end
T2=norm([x(end)-x(1);y(end)-y(2);0]);
inc_v=[x(end)-x(1);y(end)-y(2);0]/norm([x(end)-x(1);y(end)-y(2);0]);
r_v=[x(1)-x1;y(1)-y1;0]/norm([x(1)-x1;y(1)-y1;0]);
rot_dir1=-sum(sign(cross(r_v,inc_v)));

%plot 1st circles
% Alp1=ph1-pi/2;
if Alp1<0
    Alp1=Alp1+2*pi;
end
if Alp1>=B1 && rot_dir1==1 && p1==1
    tt=linspace(B1,Alp1,400);
elseif Alp1<B1 && rot_dir1==1 && p1==1
    tt=linspace(Alp1,Alp1-(2*pi-(B1-Alp1)),400);
elseif Alp1>=B1 && rot_dir1==-1 && p1==-1
    tt=linspace(B1,B1-(2*pi-(Alp1-B1)),400);
elseif Alp1<B1 && rot_dir1==-1 && p1==-1
    tt=linspace(Alp1,B1,400);    
else
     L=1e15;
    disp('not feasible')
    return
end
T1=r1*(max(tt)-min(tt));
XC1=zeros(400,2);
for i=1:1:400
   XC1(i,:)=[xxc1;yyc1]+R(th)*[r1*cos(tt(i));r1*sin(tt(i))];
end
% plot 2nd circles
% Alp=ph2-pi/2;
if Alp2<0
    Alp2=Alp2+2*pi;
end
B2=B1;
if B2<0
    B2=B2+2*pi;
end
if B2>2*pi
    B2=B2-2*pi;
end

inc_v=[x(end)-x(1);y(end)-y(2);0]/norm([x(end)-x(1);y(end)-y(2);0]);
r_v=[x(end)-x2;y(end)-y2;0]/norm([x(end)-x2;y(end)-y2;0]);
rot_dir=-sum(sign(cross(r_v,inc_v)));

if Alp2>=B2 && rot_dir==1&& p2==1
    tt=linspace(B2,Alp2-2*pi,400);
elseif Alp2<B2 && rot_dir==1 && p2==1
    tt=linspace(B2,Alp2,400);
elseif Alp2>=B2 && rot_dir==-1 && p2==-1
     tt=linspace(B2,Alp2,400);
elseif Alp2<B2 &&  rot_dir==-1 && p2==-1
    tt=linspace(B2,2*pi+Alp2,400);         
else
    L=1e15;
    disp('not feasible')
    return
end
T3=r2*(max(tt)-min(tt));
XC2=zeros(400,2);
for i=1:1:400
   XC2(i,:)=[xxc2;yyc2]+R(th)*[r2*cos(tt(i));r2*sin(tt(i))];
end
L=T1+T2+T3;
end


% now draw the arrows
XA1=[x0,y0;x0+0.9*max(r1,r2)*cos(ph0),y0+0.9*max(r1,r2)*sin(ph0)];
XA2=[xf,yf;xf+0.9*max(r1,r2)*cos(phf),yf+0.9*max(r1,r2)*sin(phf)];

plot(XC1(:,1),XC1(:,2),'r',XC2(:,1),XC2(:,2),'b',XL(:,1),XL(:,2),'k',XL(1,1),XL(1,2),'k',XL(end,1),XL(end,2),'k',XA1(:,1),XA1(:,2),'r',XA2(:,1),XA2(:,2),'b',x0,y0,'rs',xf,yf,'bs','linewidth',2,'MarkerSize',10)
axis square
axis equal






