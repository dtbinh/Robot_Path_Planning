function V=get_side_ind2(xk,XY,R)
%given xk get all the adjacent positions and manvs of cell to the current
%cell

N=size(XY,1);
x=xk(1:2);
v=xk(3:4);
v=v(:);
vn=v/norm(v);
x=x(:);
V=[];

ind=sqrt(sum((repmat(x',N,1)-XY(:,1:2)).^2,2))<=R;
xks=XY(ind,1:2);
mans=XY(ind,5);
vec=xks-repmat(x',length(mans),1);
th=zeros(length(mans),1);
for i=1:1:length(mans)
th(i)=dot(vec(i,:)/sqrt(sum(vec(i,:).^2,2)),vn);
end
V=vertcat(V,[vec,mans,xks,th]);





end