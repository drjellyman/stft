function [ ssd ] = myGetDist( xloc,sloc )
    % Returns the distances between each vector in xloc and sloc
    Nx = length(xloc(1,:));
    Ns = length(sloc(1,:));
    ssd = zeros(Ns,Nx); % ssd = source to sensor distances
    for ns=1:Ns
        for nx=1:Nx
            ssd(ns,nx) = norm(xloc(:,nx)-sloc(:,ns));
        end
    end
end

