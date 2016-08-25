function F=max_ent_surf_fit(lam,s,phi)
F=zeros(size(s,1),1);
for i=1:1:size(s,1)
F(i)=(exp(sum(lam.*phi(s(i,:)))));
end
max(F);
% F=sum(F);
end