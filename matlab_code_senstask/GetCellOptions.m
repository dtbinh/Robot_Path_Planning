function [st,M]=GetCellOptions(further,opt)
% get the parameter to the right of the options
% st isthe state if the option exists or not, M is the entry just to the
% right of this entry
ind=strcmpi(further,opt);
ind=[ind,0];
% error if the option is not there or if there two entreis with same option

if max(ind)==0 || length(ind(ind==1))>1
    M=NaN;
    st=false;
    return
end

[~,loc]=max(ind);
if loc==length(further)
    M=NaN;
    st=true;
    return
end

M=further{loc+1};
st=true;

