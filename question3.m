%%% define the parameter
%gamma = 4;
clear;
clc;
time = 100;
N = 10;
sigma = 1;
B=2;
W = [1,-0.5,-1,0.5]; 
errordata = zeros(5,size(10^-B:10^(-0.1):10^B,2));
index =1;
for gamma =  10^-B:10^(-0.1):10^B %%% iterate gamma from 0.01 to 100
    L2 = zeros(1,time);
for iteration = 1:time %%%% for every gamma iterate for 100 times
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
Westimated = LikelyYX * (LikelyXX+10*(sigma^2*(gamma^2*eye(4))^(-1)))^-1;

%%% calculate the L2 distance
L2(iteration) = (sum((Westimated - W).^ 2));
end

%%%% sort the total L2 distance
L2Sort = sort(L2,'ascend');
errordata(1,index) = L2Sort(1);
errordata(2,index) = L2Sort(time/4);
errordata(3,index) = L2Sort(time/2);
errordata(4,index) = L2Sort(3*time/4);
errordata(5,index) = L2Sort(time);
index = index+1;

%%% plot the value (1) is the min, (time/4) is 25%, (time/2) is median, (time*3/4) is
%%% 75% and (time) is the maximum 
%plot(gamma,L2Sort(1),'*');%,gamma,L2Sort(time/4),'o',gamma,L2Sort(time/2),'+',gamma,L2Sort(time*3/4),'.',gamma,L2Sort(time),'d');
%hold on
end
minstr = sprintf("min * %4.3f",max(errordata(1,:)));
qauterstr = sprintf("25 * %4.3f",max(errordata(2,:)));
medianstr = sprintf("median * %4.3f",max(errordata(3,:)));
quater3str = sprintf("75 * %4.3f",max(errordata(4,:)));
maxstr = sprintf("max * %4.3f",max(errordata(5,:)));

errordata(1,:) = errordata(1,:)./max(errordata(1,:));
errordata(2,:) = errordata(2,:)./max(errordata(2,:));
errordata(3,:) = errordata(3,:)./max(errordata(3,:));
errordata(4,:) = errordata(4,:)./max(errordata(4,:));
errordata(5,:) = errordata(5,:)./max(errordata(5,:));


gamma = 10^-B:10^(-0.1):10^B;
plot(gamma,errordata(1,:),'o')
hold on
plot(gamma,errordata(2,:),'s')
plot(gamma,errordata(3,:),'*')
plot(gamma,errordata(4,:),'+')
plot(gamma,errordata(5,:),'x')
legend(minstr,qauterstr,medianstr,quater3str,maxstr)
xlabel('gamma value','FontSize',16);
ylabel('square error','FontSize',16);
txtTitle = sprintf('simga = %f', sigma);
title(txtTitle,'FontSize',16);
