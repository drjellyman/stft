function [ y ] = myMse( Wopt,W )
    % Calculate the mean squared error between two weights matrices
    
   D = abs(Wopt-W).^2;
   y = sum(D(:))/numel(Wopt);


end

