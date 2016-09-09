function colld=collision_penalty(Xpath,i,j,Ncs,k)
if Ncs==0
    colld=0;
    return
end
if Xpath{1}(1,1)==-1
    colld=0;
    return
end

colld=zeros(Ncs,1);
sigm=2;
for ss=1:1:Ncs
    P=norm(Xpath{ss}(k,1:2)-[i,j]);
    if P==0
       colld=1e10;
       break
    end
    colld(ss)=1*exp(-0.5*P^2/sigm^2);
end

colld=max(colld);

end