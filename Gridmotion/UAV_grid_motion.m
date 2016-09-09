function [Xf,V]=UAV_grid_motion(ix,jy,dir,u,SensorModel)
% get the control actions in order
% their cost and their end point
% arrange them correctly
% return the final UAV positiona and the energy involved 

% temp{4}.J=temp{4}.J*3;
% temp{5}.J=temp{5}.J*3;
% temp{8}.J=temp{8}.J*3;
% temp{9}.J=temp{9}.J*3;
% temp{3}.J=temp{3}.J/4;
% temp{2}.J=temp{2}.J/3;
% temp{6}.J=temp{6}.J/3;
%     keyboard

N=length(SensorModel.Move.Gridx0:SensorModel.Move.Griddx:SensorModel.Move.Gridxlim);

% if u==3 || u==4 || u==8
% 
% else
%     Xf=nan;
%     V=nan;
%     return 
% end
% define boundary conditions
if (ix==1 && jy==1 && dir==3 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7||u==8||u==9))
   Xf=nan;
   V=nan;
   return 
end
if (ix==1 && jy==1 && dir==4 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7||u==4||u==5))
   Xf=nan;
   V=nan;
   return 
end
if (ix==1 && jy==N && dir==3 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7||u==4||u==5))
   Xf=nan;
   V=nan;
   return 
end
if (ix==1 && jy==N && dir==2 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7||u==8||u==9))
   Xf=nan;
   V=nan;
   return 
end
%%%
if (ix==N && jy==1 && dir==1 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7||u==4||u==5))
   Xf=nan;
   V=nan;
   return 
end
if (ix==N && jy==1 && dir==4 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7||u==8||u==9))
   Xf=nan;
   V=nan;
   return 
end
%%%
if (ix==N && jy==N && dir==1 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7||u==8||u==9))
   Xf=nan;
   V=nan;
   return 
end
if (ix==N && jy==N && dir==2 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7||u==4||u==5))
   Xf=nan;
   V=nan;
   return 
end
%%%
% on boundary pointing out
if (ix==1 && dir==3 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7))
   Xf=nan;
   V=nan;
   return
end
% on boundary pointing out
if (ix==N && dir==1 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7))
   Xf=nan;
   V=nan;
   return
end

% on boundary pointing out
if (jy==1 && dir==4 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7))
   Xf=nan;
   V=nan;
   return
end
% on boundary pointing out
if (jy==N && dir==2 && (u==1 ||u==2 ||u==3 ||u==6 ||u==7))
   Xf=nan;
   V=nan;
   return
end

% on south(4) boundary pointing east(1)
   if jy==1 && dir==1 && (u==1  || u==2 || u==4  || u==5)
   Xf=nan;
   V=nan;
   return
   end 

% on south(4) boundary and  pointing west
   if jy==1 && dir==3 && (u==7  || u==6 || u==8  || u==9)
   Xf=nan;
   V=nan;
   return
   end 
   
% on north boundary and  pointing east
   if jy==N && dir==1 && (u==7  || u==6 || u==8  || u==9)
   Xf=nan;
   V=nan;
   return
   end 
   
   % on north boundary and  pointing west
   if jy==N && dir==3 && (u==1  || u==2 || u==4  || u==5)
   Xf=nan;
   V=nan;
   return
   end
   
   % on west boundary and  pointing north
   if ix==1 && dir==2 && (u==7  || u==6 || u==8  || u==9)
   Xf=nan;
   V=nan;
   return
   end
    % on west boundary and  pointing south
   if ix==1 && dir==4 && (u==1  || u==2 || u==4  || u==5)
   Xf=nan;
   V=nan;
   return
   end
   % on east boundary and  pointing north
   if ix==N && dir==2 && (u==1  || u==2 || u==4  || u==5)
   Xf=nan;
   V=nan;
   return
   end
    % on east boundary and  pointing south
   if ix==N && dir==4 && (u==7  || u==6 || u==8  || u==9)
   Xf=nan;
   V=nan;
   return
   end

%% if non of the cases then these are most generalized case
%%%%%%%%%%%%%%%%%%%%%%%%% dir =2 %%%%%%%%%%%%%%%%%%%%%%%%%%%
if u==1 && dir==2 
   Xf(1)=ix+1;
   Xf(2)=jy+1;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==2 && dir==2
   Xf(1)=ix+1;
   Xf(2)=jy+1;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==3 && dir==2
   Xf(1)=ix;
   Xf(2)=jy+1;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==4 && dir==2
   Xf(1)=ix+1;
   Xf(2)=jy;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==5 && dir==2
   Xf(1)=ix+1;
   Xf(2)=jy;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==6 && dir==2 
   Xf(1)=ix-1;
   Xf(2)=jy+1;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
if u==7 && dir==2
   Xf(1)=ix-1;
   Xf(2)=jy+1;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==8 && dir==2
   Xf(1)=ix-1;
   Xf(2)=jy;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
if u==9 && dir==2
   Xf(1)=ix-1;
   Xf(2)=jy;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
%%%%%%%%%%%%%%%%%% dir =1 %%%%%%%%%%%%%%%%%%
if u==1 && dir==1 
   Xf(1)=ix+1;
   Xf(2)=jy-1;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==2 && dir==1
   Xf(1)=ix+1;
   Xf(2)=jy-1;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==3 && dir==1
   Xf(1)=ix+1;
   Xf(2)=jy;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==4 && dir==1
   Xf(1)=ix;
   Xf(2)=jy-1;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==5 && dir==1
   Xf(1)=ix;
   Xf(2)=jy-1;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==6 && dir==1 
   Xf(1)=ix+1;
   Xf(2)=jy+1;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
if u==7 && dir==1
   Xf(1)=ix+1;
   Xf(2)=jy+1;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==8 && dir==1
   Xf(1)=ix;
   Xf(2)=jy+1;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
if u==9 && dir==1
   Xf(1)=ix;
   Xf(2)=jy+1;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
%%%%%%%%%%%%%%%%%% dir =3 %%%%%%%%%%%%%%%%%%
if u==1 && dir==3 
   Xf(1)=ix-1;
   Xf(2)=jy+1;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==2 && dir==3
   Xf(1)=ix-1;
   Xf(2)=jy+1;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==3 && dir==3
   Xf(1)=ix-1;
   Xf(2)=jy;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==4 && dir==3
   Xf(1)=ix;
   Xf(2)=jy+1;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==5 && dir==3
   Xf(1)=ix;
   Xf(2)=jy+1;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==6 && dir==3 
   Xf(1)=ix-1;
   Xf(2)=jy-1;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
if u==7 && dir==3
   Xf(1)=ix-1;
   Xf(2)=jy-1;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==8 && dir==3
   Xf(1)=ix;
   Xf(2)=jy-1;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
if u==9 && dir==3
   Xf(1)=ix;
   Xf(2)=jy-1;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
%%%%%%%%%%%%%%%%%% dir =4 %%%%%%%%%%%%%%%%%%
if u==1 && dir==4 
   Xf(1)=ix-1;
   Xf(2)=jy-1;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==2 && dir==4
   Xf(1)=ix-1;
   Xf(2)=jy-1;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==3 && dir==4
   Xf(1)=ix;
   Xf(2)=jy-1;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==4 && dir==4
   Xf(1)=ix-1;
   Xf(2)=jy;
   Xf(3)=3;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==5 && dir==4
   Xf(1)=ix-1;
   Xf(2)=jy;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==6 && dir==4 
   Xf(1)=ix+1;
   Xf(2)=jy-1;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
if u==7 && dir==4
   Xf(1)=ix+1;
   Xf(2)=jy-1;
   Xf(3)=4;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
return
end
if u==8 && dir==4
   Xf(1)=ix+1;
   Xf(2)=jy;
   Xf(3)=1;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
if u==9 && dir==4
   Xf(1)=ix+1;
   Xf(2)=jy;
   Xf(3)=2;
   V=SensorModel.Move.GridTemplate{u}.J; % control cost
   return
end
