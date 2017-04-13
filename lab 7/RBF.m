function f=RBF(X,func,sigma)

%i=length(H);
switch func
    case 1
        for i=1:length(X)
                f(i)=exp((-X*X)/2*(sigma*sigma))
        end
    
    case 2
        for i=1:length(X)
            f(i)=X*X*log(X)
        end
    otherwise
        f(i)=0;
end

               
        
       
        
    
