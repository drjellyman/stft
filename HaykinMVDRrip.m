clear all; 
close all; 

p = 8; % no of sensors
Ninit = p; 
Nsnaps = 300; % No of iterations
Ndata = Ninit + Nsnaps;  % Total data
mean_v = 0; 
var_v = 1; 
TNRdB = 10; % Target to noise ratio
INRdB = 40; % interferer to noise ratio
numst = 301;  % spatial resolution
sin_theta = [-0.2, 0]; % T and I location
phi = pi.*sin_theta; % electrical angle of T and I
A = sqrt(var_v)*10.^([TNRdB, INRdB]./20); % weighting for T and I
e = exp(-j*[1:(p-1)]'*phi(1)); % steering vector in look direction

sig_x = A(1)*exp(j*[1:p]*phi(1)); % input signal weighting and direction by sensor
for i=1:Ndata,
    v_tmp = sqrt(var_v/2) * randn(2,p)+mean_v; % noise temp
    v = v_tmp(1,:)+j*v_tmp(2,:);    
    Psi = 2*pi*rand; 
    Xi(i,:) = sig_x + A(2)*exp(j*[1:p]*phi(2)+Psi)+v; % observation  
end

g = 1; % unity gain
d = g*Xi;   % desired signal
u = diag(Xi(:,1))*(ones(Ndata,1)*e.')-Xi(:,2:p); % 
mu = 1e-10;

% [W,xp] = lms(u,d,mu);
N = min(size(u,1), size(d,1));
Nin = size(u,2);
Nout = size(d,2);
w = zeros(Nout,Nin);
W = [];
for n = 1:N
    W = [W;w];
    xp(n,:) = u(n,:)*w'; % prediction of next sampel
    ee(n,:) = d(n,:)-xp(n,:); % prediction of error
    w = w + mu*ee(n,:)' * u(n,:); % adapt weight matrix
end

W0 = g-W*conj(e);
W = [W0,W];

WH = conj(W(Ndata,:)); % use last weight vector
st = linspace(-1,1,numst); % set up sinetheta space
est = exp(-j*pi*[0:(p-1)]'*st); % steering matrix
P = 20*log10(abs(WH*est).^2); % compute response to differennt steering vectors

R = (Xi'*Xi)/Ndata; % Correlation matrix from generated data
R_inv = inv(R);

for i = 1:numst
    s = exp(-j*pi*[0:(p-1)]'*st(i));
    P_MVDR(i) = abs(1/(s'*R_inv*s));
end
P_MVDR = 20*log10(P_MVDR./max(P_MVDR));

figure; plot(st,P)
