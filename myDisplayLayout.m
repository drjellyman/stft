function [  ] = myDisplayLayout( xloc,sloc )
    figure; plot3(xloc(1,:), xloc(2,:), xloc(3,:), '*'); grid on; hold on; 
    plot3(sloc(1,1), sloc(2,1), sloc(3,1), 'o'); 
    plot3(sloc(1,2:end), sloc(2,2:end), sloc(3,2:end), '^'); legend('Sensors','Target','Interferer')
    set(gca, 'fontsize', 14); title('Sensor and source location');
end

