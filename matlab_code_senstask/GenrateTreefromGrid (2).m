function Tk=GenrateTreefromGrid(XY,XYADJindex,K,x0ind)
% XY is 2D grid of points arranged as a vector
% XYindex is a vector (col) of index/order of the corresponding drip points
% in XY
% d is grid spacing: For now only works on equispaced sqaure grids on x-y
% space
% x0 is the location to start the tree from
% 1:K levels
% XYADJindex{} is a cell with adjacend node index for the given index node

% The tree is represented as :
% at level k
% Tk=cell(nNodes,3):  Tk{1,1}=18 (curr node number),   Tk{1,2}=Left
% connected nodes,  Tk{1,3}= Right connected ndoes

ng=size(XY,1);
% x0=x0(:)';

Tk=cell(K,1);
for i=1:1:K
    Tk{i}=cell(1,3);
end

% S=repmat(x0,ng,1)-XY;
% 
% [~,x0ind]=min(sqrt(sum(S.^2,2)));

Tk{1}{1,1}=x0ind;

for i=1:1:K-1
    if i==1 % then only right chldered
        for j=1:1:size(Tk{i},1)
            Tk{i}{j,2}=[];
            xind=Tk{i}{j,1};
            Tk{i}{j,3}=XYADJindex{xind};
            Tk{i+1}=UpdateRightTree(Tk{i+1},xind,Tk{i}{j,3});
            %             R=horzcat(R,Tk{i}{j,3});
        end
        %         R=unique(R);
        %         for r=1:1:length(R)
        %             Tk{i+1}{r,1}=R(r);
        %         end
        %         for j=1:1:size(Tk{i},1)
        %             xind=Tk{i}{j,1};
        %             for r=1:1:length(Tk{i}{j,3})
        %              Tk{i}{j,3}(r)
        %         end
        continue
    end

    for j=1:1:size(Tk{i},1)
        xind=Tk{i}{j,1};
        Tk{i}{j,3}=XYADJindex{xind};
        Tk{i+1}=UpdateRightTree(Tk{i+1},xind,Tk{i}{j,3});

    end
    
end

end

% first create the nodes rk at time k1
% take the currenet node xk and put into left of k1 node
function RT=UpdateRightTree(RT,xk,rk)

rk=unique(rk);
if size(RT{1,1},1)==0
    for i=1:1:length(rk)
        RT{i,1}=rk(i);
        RT{i,2}=xk;
    end
    return
end
% keyboard
for i=1:1:length(rk)
    flg=0;
    for k=1:1:size(RT,1)
        if RT{k,1}==rk(i)
            RT{k,2}=horzcat(RT{k,2},xk);
            flg=1;
        end
    end
    if flg==0  % if not there
        kk=size(RT,1)+1;
        RT{kk,1}=rk(i);
        RT{kk,2}=xk;
    end
end


end