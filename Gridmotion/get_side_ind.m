function V=get_side_ind(xk,XY,dx)
%given xk get all the adjacent positions and manvs of cell to the current
%cell

N=size(XY,1);
x=xk(1:2);
v=xk(3:4);
v=v(:);
vn=v/norm(v);
x=x(:);
V=[];

[a,bc]=min(sqrt(sum((repmat(x',N,1)-XY(:,1:2)).^2,2)));
xks=XY(bc,1:2);
xc=xks;
mans=XY(bc,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);




% right
[~,s]=min(sqrt(sum((repmat(xc+[dx,0],N,1)-XY(:,1:2)).^2,2)));
xks=XY(s,1:2);
mans=XY(s,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);

% right up
[~,s]=min(sqrt(sum((repmat(xc+[dx,dx],N,1)-XY(:,1:2)).^2,2)));
xks=XY(s,1:2);
mans=XY(s,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);

% right down
[~,s]=min(sqrt(sum((repmat(xc+[dx,-dx],N,1)-XY(:,1:2)).^2,2)));
xks=XY(s,1:2);
mans=XY(s,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);


% left
[~,s]=min(sqrt(sum((repmat(xc+[-dx,0],N,1)-XY(:,1:2)).^2,2)));
xks=XY(s,1:2);
mans=XY(s,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);


% left up
[~,s]=min(sqrt(sum((repmat(xc+[-dx,dx],N,1)-XY(:,1:2)).^2,2)));
xks=XY(s,1:2);
mans=XY(s,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);


% left down
[~,s]=min(sqrt(sum((repmat(xc+[-dx,-dx],N,1)-XY(:,1:2)).^2,2)));
xks=XY(s,1:2);
mans=XY(s,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);


% up
[~,s]=min(sqrt(sum((repmat(xc+[0,dx],N,1)-XY(:,1:2)).^2,2)));
xks=XY(s,1:2);
mans=XY(s,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);


% down
[~,s]=min(sqrt(sum((repmat(xc+[0,-dx],N,1)-XY(:,1:2)).^2,2)));
xks=XY(s,1:2);
mans=XY(s,5);
vec=xks-x';
th=dot(vec/norm(vec),vn);
V=vertcat(V,[vec,mans,xks,th]);

% ind = find(V(:,end)>=0);
% V=V(ind,:);



end