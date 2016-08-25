function plot_dubinspath(r0,rf,X0,Xf,p1,p2,dir)

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
th=atan2(y2,x2);
Xx=R(-th)*[x2;y2];
x2=Xx(1);
y2=Xx(2);
ph1=ph0-th;
ph2=phf-th;
d1=norm([x2,y2]);


if dir==1
B1=atan(sqrt(r1^2*d1^2/(r1+r2)^2-r1^2)/r1);

% first plot two circles
tt=linspace(0,2*pi,200);
XC1=zeros(200,2);
XC2=zeros(200,2);
for i=1:1:200
   XC1(i,:)=[xxc1+r1*cos(tt(i)),yyc1+r1*sin(tt(i))];
   XC2(i,:)=[xxc2+r2*cos(tt(i)),yyc2+r2*sin(tt(i))];
end

% line between the contacr points
x=(r1*cos(B1)):0.1:(d1+r2*cos(B1+pi));
y=(r2*sin(B1+pi))+(r1*sin(B1)-r2*sin(B1+pi))*(x-d1-r2*cos(B1+pi))/(r1*cos(B1)-d1-r2*cos(B1+pi));
XL=zeros(length(x),2);
for i=1:1:length(x)
   XL(i,:)=([xxc1;yyc1]+R(th)*[x(i);y(i)])'; 
end
end

if dir==2
B1=-atan(sqrt(r1^2*d1^2/(r1+r2)^2-r1^2)/r1);

% first plot two circles
tt=linspace(0,2*pi,200);
XC1=zeros(200,2);
XC2=zeros(200,2);
for i=1:1:200
   XC1(i,:)=[xxc1+r1*cos(tt(i)),yyc1+r1*sin(tt(i))];
   XC2(i,:)=[xxc2+r2*cos(tt(i)),yyc2+r2*sin(tt(i))];
end

% line between the contacr points
x=(r1*cos(B1)):0.1:(d1+r2*cos(B1+pi));
y=(r2*sin(B1+pi))+(r1*sin(B1)-r2*sin(B1+pi))*(x-d1-r2*cos(B1+pi))/(r1*cos(B1)-d1-r2*cos(B1+pi));
XL=zeros(length(x),2);
for i=1:1:length(x)
   XL(i,:)=([xxc1;yyc1]+R(th)*[x(i);y(i)])'; 
end
end

if dir==3
B1=atan(sqrt(d1^2-(r1-r2)^2)/abs(r1-r2));

% first plot two circles
tt=linspace(0,2*pi,200);
XC1=zeros(200,2);
XC2=zeros(200,2);
for i=1:1:200
   XC1(i,:)=[xxc1+r1*cos(tt(i)),yyc1+r1*sin(tt(i))];
   XC2(i,:)=[xxc2+r2*cos(tt(i)),yyc2+r2*sin(tt(i))];
end

% line between the contacr points
x=(r1*cos(B1)):0.1:(d1+r2*cos(B1));
y=(r2*sin(B1))+(r1*sin(B1)-r2*sin(B1))*(x-d1-r2*cos(B1))/(r1*cos(B1)-d1-r2*cos(B1));
XL=zeros(length(x),2);
for i=1:1:length(x)
   XL(i,:)=([xxc1;yyc1]+R(th)*[x(i);y(i)])'; 
end
end

if dir==4
B1=-atan(sqrt(d1^2-(r1-r2)^2)/abs(r1-r2));

% first plot two circles
tt=linspace(0,2*pi,200);
XC1=zeros(200,2);
XC2=zeros(200,2);
for i=1:1:200
   XC1(i,:)=[xxc1+r1*cos(tt(i)),yyc1+r1*sin(tt(i))];
   XC2(i,:)=[xxc2+r2*cos(tt(i)),yyc2+r2*sin(tt(i))];
end

% line between the contacr points
x=(r1*cos(B1)):0.1:(d1+r2*cos(B1));
y=(r2*sin(B1))+(r1*sin(B1)-r2*sin(B1))*(x-d1-r2*cos(B1))/(r1*cos(B1)-d1-r2*cos(B1));
XL=zeros(length(x),2);
for i=1:1:length(x)
   XL(i,:)=([xxc1;yyc1]+R(th)*[x(i);y(i)])'; 
end
end


% now draw the arrows
XA1=[x0,y0;x0+0.9*max(r1,r2)*cos(ph0),y0+0.9*max(r1,r2)*sin(ph0)];
XA2=[xf,yf;xf+0.9*max(r1,r2)*cos(phf),yf+0.9*max(r1,r2)*sin(phf)];

plot(XC1(:,1),XC1(:,2),'r',XC2(:,1),XC2(:,2),'b',XL(:,1),XL(:,2),'k',XL(1,1),XL(1,2),'ko',XL(end,1),XL(end,2),'ko',XA1(:,1),XA1(:,2),'r',XA2(:,1),XA2(:,2),'b',x0,y0,'rs',xf,yf,'bs','linewidth',2,'MarkerSize',10)
axis square
axis equal







end