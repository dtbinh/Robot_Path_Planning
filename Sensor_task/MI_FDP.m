function MI=MI_FDP(P0,PKU,flag)

MI=0;
for ntar=1:1:length(P0)
   

                if det(P0{ntar})/det(PKU{ntar})<1.01  %Pf/Pu
                    MI=MI+0;
                elseif isnan(flag(ntar))
                    MI=MI+0;
                else
                    MI=MI+0.5*log(det(P0{ntar})/det(PKU{ntar}));
                end

    
end





end