function [X,Extrap2ptErr] = twoPtExtrap(dx,values,exact)
%Preforms 3-point extrapolation and produces error values 
%   vector of dx values (Delta X)
%   vector of values corresponding to dx values for a given postion
%   exact solution for given position


   for m = 1:length(dx)-1
       
       X(m) = dx(m);
            
       Ev2point(m) = (4*values(m+1) - values(m))/3;
       
       Extrap2ptErr(m) = abs(exact - Ev2point(m))./(exact*100);
   end






end