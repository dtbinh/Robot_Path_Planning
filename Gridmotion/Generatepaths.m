function P=Generatepaths(Tk,p,k)
% p is parent node at k-1

P=[];

for j=1:1:size(Tk{k},1)
    if Tk{k}{j,1}==p
        if isempty(Tk{k}{j,3})
            P=p;
            return
        end
        
        for h=1:1:length(Tk{k}{j,3})     
            S=Generatepaths(Tk,Tk{k}{j,3}(h),k+1);
            P=vertcat(P,[p*ones(size(S,1),1),S]);
        end
        break
    end
end
















end

