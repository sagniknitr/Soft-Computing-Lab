clc;
clear all;
close all;
h1=[1 0.5];
it=1:100;
N=1000;
mse=zeros(1,N); 
e=zeros(N,length(it));
e1=e;
for k=1:length(it)
x11=randsrc(1,N+1);
y11=zeros(1,length(x11));
y12=y11;
x12=y11; 
for i=2:length(x11)
xin=[x11(i) x11(i-1)]; 
y11(i)=h1*xin'; 
y12(i)=awgn(y11(i),15);
end

w11=rand(1,2);
b=rand(1,1); 
for i=2:length(y12)
yin=[y12(i) y12(i)];% d=0 delay elementsequalizer 
x12(i)=hardlims(w11*yin'+b); 
e(i-1,k)=x11(i)-x12(i); 
w11=w11+0.05*e(i-1,k).*yin;
b=b+0.05*e(i-1,k);
end
end
eit=zeros(1,length(y12)-1); 
for i=2:length(y12)
for k=1:length(it) 
eit(i-1)=eit(i-1)+(e(i-1,k)^2);
end
 eit(i-1)=eit(i-1)/(length(it));
end
figure();

plot(2:length(x11),y12(2:length(x11))),grid on,title('Observed Channel Output'); xlabel('Input'),ylabel('Channel Output');
figure();
semilogy(1:length(eit),eit),grid on,title('MSE vs iterations'),xlabel('Iterations'),ylabel('MSE'); figure();
plot([-2:2],-1*(w11(1)*[-2:2]+b)/w11(2)),grid on,title('Decision Boundary'); xlabel('X(n)'),ylabel('X(n-1)'),hold on
for i=2:length(x11) 
    xin=[x11(i) x11(i-1)]; 
    if(x11(i)==1)
        plot(x11(i),x11(i-1),'*r');
    else
        plot(x11(i),x11(i-1),'or');
    end
y12=awgn(xin,20); 
if(x11(i)==1)
    plot(y12(1),y12(2),'*b'); 
else
plot(y12(1),y12(2),'ob'); 
end
end