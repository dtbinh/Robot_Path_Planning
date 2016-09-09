function [SensorModel,Trajectories]=UAV_traj_pplied_control(ss,SensorModel,Trajectories,curk)
% from curk to curk+1
% position at curk+1 is updated

c0=SensorModel.Move.Gridpos{ss}(curk,:);
u=SensorModel.Move.DP_MU{ss,curk}(c0(1),c0(2),c0(3));
SensorModel.Move.GridOptiControl{ss}(curk)=u;

% [Xf,V]=UAV_grid_motion(c0(1,1),c0(1,2),c0(1,3),u,xlim,d,templates1);
% try
if u==0 || isempty(u) || isnan(u)
    disp('incompatible u')
    keyboard
end
[Xf,~]=SensorModel.Move.GridMotion(c0(1,1),c0(1,2),c0(1,3),u);
% catch
%     keyboard
% end
XX=SensorModel.Move.GridTemplate2{u,1}(:,[1,2]);
phi=SensorModel.Move.GridTemplate2{u,1}(:,3);
            Y=XX;
            if c0(1,3)==1  %right
                Y(:,1)=XX(:,2);
                Y(:,2)=-XX(:,1);
                phi=phi-pi/2;
            end
            if c0(1,3)==2  %up , the default
                Y=XX;
            end
            if c0(1,3)==3 %left
                Y(:,1)=-XX(:,2);
                Y(:,2)=XX(:,1);
                phi=phi+pi/2;
            end
            if c0(1,3)==4   %down
                Y(:,1)=-XX(:,1);
                Y(:,2)=-XX(:,2);
                phi=phi+pi;
            end
%             keyboard
            Trajectories{ss}=vertcat(Trajectories{ss},[Y+repmat([SensorModel.Move.Xgrid(c0(1,1),c0(1,2)),SensorModel.Move.Ygrid(c0(1,1),c0(1,2))],size(SensorModel.Move.GridTemplate2{u,1},1),1),phi,curk*ones(length(phi),1)]);
            %           keyboard
            c0(1,1)=Xf(1);c0(1,2)=Xf(2);c0(1,3)=Xf(3);
            SensorModel.Move.Gridpos{ss}(curk+1,:)=c0;
end