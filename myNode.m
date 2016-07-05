classdef myNode
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        loc         % 3d location of node
        W           % Local weights matrix
        L           % Local Lambda dual weight matrix
        R           % Observation covariance
        d           % Look direction vector
        N           % Neighbor ids
        Nlen        % Number of neighbors including this node and it's virtual node
        Amn         % Consistency matrix from this node (m) to node n
%         Anm         % Consistency matrix from node n to node m
        X           % Observation data in frequency domain for the neighborhood of this node
        Y           % MVDR output of node m
    end
    
    methods
    end
    
end

