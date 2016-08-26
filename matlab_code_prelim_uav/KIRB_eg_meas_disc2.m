 function yk=KIRB_eg_meas_disc2(x,para)
%     zi=x(1)-0;
%     ni=x(3)+0;
   
%     yk(1,1)=sqrt(x(1)^2+x(2)^2);
    yk(2,1)=atan2(x(2),x(1));
yk(1,1)=sqrt((x(1))^2+(x(2))^2);
%  yk(2,1)=sqrt((x(1)-1e4)^2+(x(2)+1.6e4)^2);
%  yk(2,1)=para(1);
%     yk(1,1)=x(1);
%     yk(2,1)=x(2);
    end