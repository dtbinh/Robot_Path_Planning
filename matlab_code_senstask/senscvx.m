function f=senscvx(x,Px,Pxz,Pz)
v1=x(1);
v2=x(2);
v3=x(3);

PPxz=Pxz;
    PPxz(:,1)=PPxz(:,1)*v1;
    PPxz(:,2)=PPxz(:,2)*v2;
    PPxz(:,3)=PPxz(:,3)*v3;
    
    PPz=Pz;
    PPz(1,2)=PPz(1,2)*v1;
    PPz(1,3)=PPz(1,3)*v1;
    PPz(2,1)=PPz(2,1)*v1;
    PPz(3,1)=PPz(3,1)*v1;
    
    PPz(1,2)=PPz(1,2)*v2;
    PPz(2,1)=PPz(2,1)*v2;
    PPz(2,3)=PPz(2,3)*v2;
    PPz(3,2)=PPz(3,2)*v2;
    
    PPz(1,3)=PPz(1,3)*v3;
    PPz(3,1)=PPz(3,1)*v3;
    PPz(2,3)=PPz(2,3)*v3;
    PPz(3,2)=PPz(3,2)*v3;
Info=Px-PPxz*inv(PPz)*PPxz';
    f=trace(Info);


end