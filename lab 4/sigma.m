function f=sigma(X)

%i=length(H);
for i=1:length(X)
    f(i)=1/(1+exp(-X(i)));
end

    