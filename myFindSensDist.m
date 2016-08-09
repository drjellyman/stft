function [ sensd ] = myFindSensDist( xloc )
    M = length(xloc(1,:));
    sensd = zeros(M,M); % sensd = sensor to sensor distance
    for m=1:M
        for mm=m:M
            sensd(m,mm) = norm(xloc(:,m)-xloc(:,mm));
        end
    end
    sensd = sensd + triu(sensd,1).'; % Convert from upper triangular to symmetric

end

