function [mk1,Pk1]=KalmanUpdate_modf_angle(mk,Pk,mz,Pz,Pcc,ym,ytp)
             K=Pcc/Pz;
%              keyboard
             if max(isnan(ym))==1 || max(isempty(ym))==1
                 mk1=NaN;
                 error('no measurement here')
             else
                 d=ym(:)-mz(:);
                 ind=find(ytp==1);
                 qz=angle2quat(mz(ind),zeros(length(ind),1),zeros(length(ind),1));
                 qm=angle2quat(ym(ind),zeros(length(ind),1),zeros(length(ind),1));
                 d(ind)=quat2angle(quatmultiply(qm,quatinv(qz)));
%                  keyboard
               mk1=mk+K*(d);
%              mk1=mk+K*(ym(:)-mz(:));
%              [mk1,mk,K*(ym(:)-mz(:))]
             end
             %update
             Pk1=Pk-K*Pz*K';

end
