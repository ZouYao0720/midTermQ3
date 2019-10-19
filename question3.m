%%% define the parameter
%gamma = 4;
clear;
clc;
time = 100;
N = 10;
sigma = 10;
B=2;
W = [1,-0.5,-1,0.5]; 
for gamma = 10^-B:10^(-0.1):10^B %%% iterate gamma from 0.01 to 100
for iteration = 1:time %%%% for every gamma iterate for 100 times
    L2 = zeros(1,time);
    x = -1 + 2*rand(N,1);
    X = [x.*x.*x,x.*x,x,ones(N,1)];
    Y = X*W' + sigma*randn(N,1);
    LikelyXX = zeros(4,4);
    LikelyYX = zeros(1,4);
%%% form the likelyhood matrix
for n = 1:N
    XX(:,:,n) = X(n,:)' *X(n,:);
    LikelyXX = LikelyXX+XX(:,:,n);
    YX(n,:) = Y(n).*X(n,:);
    LikelyYX = LikelyYX + YX(n,:);
end
%%% doing the estiamting process
%%% . formular W(map) = sum(y*x) * invers(sigmg1^2*inver(cov(W)) + sum(x*x'))
Westimated = LikelyYX * ((sigma^2*(gamma^2*eye(4))^(-1))+LikelyXX)^-1;

%%% calculate the L2 distance
L2(iteration) = (sum((Westimated - W) .^ 2));
end

%%%% sort the total L2 distance
L2Sort = sort(L2,'ascend');

%%% plot the value (1) is the min, (time/4) is 25%, (time/2) is median, (time*3/4) is
%%% 75% and (time) is the maximum 
plot(gamma,L2Sort(1),'*',gamma,L2Sort(time/4),'o',gamma,L2Sort(time/2),'+',gamma,L2Sort(time*3/4),'.',gamma,L2Sort(time),'d');
hold on
end
legend('the min','25%','median','75%','max')
xlabel('gamma value','FontSize',16);
ylabel('square error','FontSize',16);
