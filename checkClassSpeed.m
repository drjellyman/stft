close all; clear

results = zeros(4,1);

%% Set up time trial for node property access
node = myNode;
node.W = rand(1000,1000);

tic
for a = 1:1000
    x = node.W(:,a);
end
results(1) = toc;

%% Set up time trial for matrix access
W = rand(1000,1000);

tic
for a = 1:1000
    x = W(:,a);
end
results(2) = toc;

%% Time trial for node property set
node = myNode;

tic
for a = 1:1000
    node.W(:,1) = rand(1000,1);
end
results(3) = toc;

%% Time trial for matrix set
tic 
for a = 1:1000
    W(:,a) = rand(1000,1);
end
results(4) = toc 
difference = [diff(results(1:2));diff(results(3:4))]



