% teting exhaustive searching method

dx=10;
B=100;
[X,Y]=meshgrid(0:dx:B);
XY=[Y(:),X(:)];
XYindex=[1:1:size(XY,1)]';
% [XY,XYindex]
ng=length(XYindex);
XYADJindex=cell(ng,1);

for i=1:1:length(XYindex)
    xy=XY(i,:);
    p=[];
    
    r=xy+[dx,0];
    if r(1)>B || r(2) >B || r(1)<0 || r(2)<0
        
    else
        [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
        p=horzcat(p,ind(1));
    end
    
    r=xy+[0,dx];
    if r(1)>B || r(2) >B || r(1)<0 || r(2)<0
        
    else
        [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
        p=horzcat(p,ind(1));
    end
    
    r=xy+[-dx,0];
    if r(1)>B || r(2) >B || r(1)<0 || r(2)<0
        
    else
        [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
        p=horzcat(p,ind(1));
    end
    
    r=xy+[0,-dx];
    if r(1)>B || r(2) >B || r(1)<0 || r(2)<0
        
    else
        [~,ind]=min(sqrt(sum((repmat(r,ng,1)-XY).^2,2)));
        p=horzcat(p,ind(1));
    end
    
    XYADJindex{i}=sort(p);
end

Tk=GenrateTreefromGrid(XY,XYADJindex,6,[50,50]);

P=Generatepaths(Tk,61,1);

S=[74,63,52]';
Pp=PrunePaths(P,S);


for i=1:1:size(Pp,1)
figure(1)

plot(XY(:,1),XY(:,2),'bo')
hold on
plot(XY(S,1),XY(S,2),'ks')

np=Pp(i,:);
plot(XY(np,1),XY(np,2),'r','linewidth',2)
pause
clf
end




