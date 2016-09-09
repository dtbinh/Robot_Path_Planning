% temp{1}=templates1{1};
% temp{2}=templates1{2};
% temp{3}=templates1{11};
% temp{4}=templates1{6};
% temp{5}=templates1{7};
% temp{6}=templates1{2};
% temp{7}=templates1{9};
% temp{8}=templates1{16};
% temp{9}=templates1{7};

[Ygrid,Xgrid]=meshgrid(0:40:400);
Vt=0;
Xf=[1,1,2];
figure(1)
while 1


if isnan(Xf)==1
    Xf=Xp;
    disp('control not allowed')
else
if Xf(3)==1
    plot(Xgrid(Xf(1),Xf(2)),Ygrid(Xf(1),Xf(2)),'ko',[Xgrid(Xf(1),Xf(2)),Xgrid(Xf(1),Xf(2))+10],[Ygrid(Xf(1),Xf(2)),Ygrid(Xf(1),Xf(2))],'b','linewidth',2);
elseif Xf(3)==2
    plot(Xgrid(Xf(1),Xf(2)),Ygrid(Xf(1),Xf(2)),'ko',[Xgrid(Xf(1),Xf(2)),Xgrid(Xf(1),Xf(2))],[Ygrid(Xf(1),Xf(2)),Ygrid(Xf(1),Xf(2))+10],'b','linewidth',2);
elseif Xf(3)==3
    plot(Xgrid(Xf(1),Xf(2)),Ygrid(Xf(1),Xf(2)),'ko',[Xgrid(Xf(1),Xf(2)),Xgrid(Xf(1),Xf(2))-10],[Ygrid(Xf(1),Xf(2)),Ygrid(Xf(1),Xf(2))],'b','linewidth',2);
elseif Xf(3)==4
    plot(Xgrid(Xf(1),Xf(2)),Ygrid(Xf(1),Xf(2)),'ko',[Xgrid(Xf(1),Xf(2)),Xgrid(Xf(1),Xf(2))],[Ygrid(Xf(1),Xf(2)),Ygrid(Xf(1),Xf(2))-10],'b','linewidth',2);
end
  axis([-10,410,-10,410]) 
  grid on
  title(num2str(Vt))
Xp=Xf;
end
Xf
  u=input('give control command');
[Xf,V]=UAV_grid_motion(Xf(1),Xf(2),Xf(3),u,400,40,templates1);
if isnan(Xf)~=1
Vt=Vt+V;
end
end
    
    
    