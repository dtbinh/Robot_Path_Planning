function Pp=PrunePaths(P,S)
% S is the set of nodes that have to be there
S=unique(S);

II=ones(size(P));
FFInd=[];
for i=1:1:length(S)
    
    FFInd=vertcat(FFInd,find(min(abs(P-S(i)*II),[],2)==0));
    FFInd=unique(FFInd);
    
end


Pp=P(FFInd,:);


