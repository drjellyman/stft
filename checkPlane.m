close all; clear all;
ns = 20;
zPos = [[1:0.3:1+(ns-1)*0.3];ones(1,ns)];
sPos = [1,10]';
for a = 1:ns
    d(a) = norm(zPos(:,a)-sPos);
    
end
figure; plot(zPos(1,:),zPos(2,:),'o',sPos(1),sPos(2),'*'); xlim([0,7]);ylim([0,11]);legend('Sensors','Source'); set(gca,'fontsize',16)
figure; plot(d); xlabel('Sensor number'); ylabel('Distance to source (m)'); set(gca,'fontsize',16)