function f=xor_test(X)

%i=length(H);

f=zeros(X,4);
for i=1:length(X)
    for j=1:4
     f(mod(j,4))=1;
     %f(ceil(3*rand(1,1)))=1;
     f(4)=xor(f(3),xor(f(1),f(2)));
    end
    
end