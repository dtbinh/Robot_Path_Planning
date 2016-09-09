function J=SensorOverlapCost(s,Xpath,SensorModel)
%s is sensor number
% Xpath is trajectory of the previous sensots (remember that these are grid coordinates)
%all the gridpoints that are being overlapped and the numbe of times
keyboard

grds=[SensorModel.Move.Xgrid(:),SensorModel.Move.Ygrid(:)];


Statgridview=cell(SensorModel.Static.Ns,1);
for i=1:1:SensorModel.Static.Ns
    xc=SensorModel.Static.pos{i}(1)
    yc=SensorModel.Static.pos{i}(2)
    
    alphaa=SensorModel.Static.FOV{i}(1)
    a=SensorModel.Static.FOV{i}(2)
   th=SensorModel.Static.FOV{i}(3)
    th=th-pi/2;
    R=[cos(th),-sin(th);sin(th),cos(th)];
    x=[0,a*tan(alphaa),-a*tan(alphaa),0];
    y=[0,a,a,0];
    % keyboard
    for j=1:1:length(x)
        xy=R*[x(j);y(j)]+[xc;yc];
        x(j)=xy(1);
        y(j)=xy(2);
    end
    
    [IN,ON] = inpolygon(grds(:,1),grds(:,2),x,y);
    
    ind=[IN;ON];
    Statgridview{i}=grds(ind==1,:);
end

% Now for moving sensors




end


function PolyonCircle(xc,yc,R,n)
% n is the number of vertices of the polygon required
th=linspace(0,2*pi,n);

end

