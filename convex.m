%%

x1 = [0:1e-3:2].';
x2 = x1;
X = x1*x2.';
y = x1./x2;
figure; surf(X,y)