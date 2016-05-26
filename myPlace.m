function [ zPos, sPos ] = myPlace( NSources, M, ag)

    % Places sensors and sources in a 3x5x8m room and returns location
    % vectors. ag = array geometry
    
    if ag == 'linear'
        zPos = repmat([1,1.5,0.5]',1,M);
        zPos(1,:) = rand(1,M)*M;
        
    end
    
    sPos = rand(3,NSources).*repmat([3,5,8]',1,NSources);
    
end

