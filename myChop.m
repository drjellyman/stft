function [ M ] = myChop( M,m)
    % Remove -1 values from M and N. Add each node to its own vector. 
    M = [M(1:min(find(M==-1))-1);m];    
end

