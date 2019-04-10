function [X,Extrap3ptErr] = threePtExtrap(dx,values,exact)
%Preforms 3-point extrapolation and produces error values 
%   vector of dx values (Delta X)
%   vector of values corresponding to dx values for a given postion
%   exact solution for given position


   for m = 1:length(dx)-3
       
       X(m) = dx(m);
      
       Ev(m) = ((values(m)*values(m+2)) - (values(m+1))^2)./(values(m)+values(m+2) - (2*values(m+1)));
       
       Extrap3ptErr(m) = abs((exact - Ev(m)))./(exact*100);

   end
   






end

