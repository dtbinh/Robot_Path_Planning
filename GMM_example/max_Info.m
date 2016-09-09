function ot=max_Info(xs,GMM,model)
% function to maximize the sensor paramters
model.hx=@(x,para)model.hx(x,para,xs);
model.gx=@(x,para)model.gx(x,para,xs);

ang=[-pi+1*pi/180:1*pi/180:pi,pi/2];
% alp=[5*pi/180,10*pi/180,20*pi/180];
alp=30*pi/180;
DKL_gsmodp=0;
for phi1=ang
     for alpha1=alp
        DKL_gsmod=ExptKL_modf_MeasOPT_gsm(GMM,model,[phi1,alpha1]);
        if DKL_gsmod>DKL_gsmodp && DKL_gsmod~=Inf
        theta=[phi1,alpha1];
        DKL_gsmodp=DKL_gsmod;
        end
     end
end
ot=[DKL_gsmodp,theta];
%  keyboard
end