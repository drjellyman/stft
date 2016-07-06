function [ Amn] = myAconsistency( M,n)

    Ml = length(M);

    Amn = [1,zeros(1,Ml-1);(M==n)];


end

