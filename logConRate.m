function [conRate] = logConRate(X,Y)
%Takes log of x and y, returns convergence rate

 func = polyfit(log(X),log(Y),1);
 
 conRate = func(1); 

end

