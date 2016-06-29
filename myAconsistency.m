function [ Amn,Anm ] = myAconsistency( M,N)
% Takes two neighbor address vectors M and N, where M is the current node 
% and N is the neighbor,and creates two A matrices that
% enforce consistency amongst the local constraints. Both A matrices are 
% 2 x k where k is the number of neighbors

    Ml = length(M);
    Nl = length(N);

    % Create Amn by stepping through N
    AmnTmp = zeros(1,Ml);
    for Ni = 1:Nl
        AmnTmp = AmnTmp + (N(Ni)==M).';
    end
    Amn = [AmnTmp;-AmnTmp];

    % Create Anm by stepping through M
    AnmTmp = zeros(1,Nl);
    for Mi = 1:Ml
        AnmTmp = AnmTmp + (M(Mi)==N).';
    end
    Anm = [-AnmTmp;AnmTmp];
end

