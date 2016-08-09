% biadmm1class7a.m is a start from basics approach to getting the perfect
% covariance working for one source and two sensors.

clear; close all;
Ns = 2; 
sLoc = [1,1,1].';
xLoc = [1,1.5,1 ; 1,1.5,1.5 ].';
t = [0:0.00001:0.005].';
f = 2e3; % (Hz)
s = exp(-1i*2*pi*f*t);
c = 343; % (m/s)
R = zeros(Ns);
sxDist = [norm(xLoc(:,1)-sLoc) ; norm(xLoc(:,2)-sLoc)];
x1 = s .* exp(-1i*2*pi*sxDist(1)/c);
x2 = s .* exp(-1i*2*pi*sxDist(2)/c);
X = [x1,x2];
for a=1:Ns
    for b=1:Ns
        R(a,b) = exp(-1i*2*pi*f*((sxDist(a)-sxDist(b)))/c) ;
    end
end
for a = 1:Ns
    d(a) = exp(-1i*2*pi*f*sxDist(a)/c);
end
d = d.';
R = R+1e-3*eye(Ns);
Wopt = (R\d)/(d'/R*d);

figure; plot(real(s)); hold on; grid on; plot(real(x1)); plot(real(x2)); legend('s','x1','x2')

y = Wopt(1)'*x1 + Wopt(2)'*x2;
figure; plot(real(y)); hold on ; plot(real(s)); grid on; legend('y','s');