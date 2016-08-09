function [  ] = myDisplayLayout( xloc,sloc )
%     dx = 0.05;
%     text = num2str([1:2]);
%     textcell = cellstr(text);
    
    figure; plot3(xloc(1,:), xloc(2,:), xloc(3,:), '*'); grid on; hold on; 
    plot3(sloc(1,1), sloc(2,1), sloc(3,1), 'o'); 
    plot3(sloc(1,2:end), sloc(2,2:end), sloc(3,2:end), '^'); legend('Sensors','Target','Interferer')
%     text(sloc(1,2:end), sloc(2,2:end), sloc(3,2:end), textcell)
    set(gca, 'fontsize', 14); title('Sensor and source location');
end

