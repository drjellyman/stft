% f*(y)=supx y^Tx-fx

x = [-2:0.01:2].';
y = x;


for yi = 1:length(y)
    fy1(yi) = max(y(yi)*x-x.^2);
end

fy2 = y.'*x-x.^2;

figure; plot(y,fy1); grid on; xlabel('y'); ylabel('f*(y)');
figure; plot(y,fy2); grid on; xlabel('y'); ylabel('f*(y)');