clc;
clear all;
close all;
input=[1 1 1 1 1 1 1 1 -1 -1 -1 -1 -1 -1 -1 -1 ; 1 1 1 1 -1 -1 -1 -1 1 1 1 1 -1 -1 -1 -1; 1 1 -1 -1 1 1 -1 -1 1 1 -1 -1 1 1 -1 -1;1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 1 -1 ];
%b=1;
%input=[1 1 1 1 -1 -1 -1 -1; 1 1 -1 -1 1 1 -1 -1;1 -1 1 -1 1 -1 1 -1];
output=[0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
weight=zeros(500,2,20);
bias = zeros(500,1,20); 
err = zeros(500,1,20);
x=input;
t=output;
%[row col]= size(x);
%out = zeros(1,4);
for k=1:20
    weight=rand([1,4])*2-1;
    bias=rand([1,1])*2-1;
    weight1(1,:,1)=weight;
    bias1(1,:,1)=bias;
    for j=1:500
 
        r = randi(4);
        x1(:,j)=x(:,r);
        y(1,j)=(weight*x(:,r)+bias);
        out(1,j) = (1/(1+exp(-y(1,j))));
        e=t(r)-out(j);
        bias=bias+e;
        weight=weight+e.*transpose(x(:,r));
        weight1(j,:,k)= weight;
        bias1(j,1,k) = bias;
        err(j,1,k) = e;
 
    end
    plot(weight1(:,1,k),bias1(:,1,k)); hold on; 
end
 
for j= 1:500
    ave = err(j,1,1);
    for k =1:20 
        ave = ave + (err(j,1,k)/20);
    end
    mse(j,1)=ave;
end
[ x1' out']
figure;
 plotpv(x,t);
 plotpc(weight,bias);
figure;
plot(mse)
