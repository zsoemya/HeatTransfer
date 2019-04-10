function [A] = KCoeffcients4thFDM(dx,b_i,k_i,nx)


%Calculate First Two Terms of Matrix, Excludes i-1 term 
beta_1 = (dx(1)^2/12)*(b_i(1)/k_i(1));

K_11 = (k_i(1)/dx(1))*(1/(1+beta_1)) + (dx(1)/2)*b_i(1); 
K_12 = (-k_i(1)/dx(1))*(1/(1+beta_1)); 

%Add to Global Matrix 
    A(1,1) = K_11; 
    A(1,2) = K_12;
    
    
%Generate all internal Global Matrix Values 
for n =2:nx

%Cycles for all values of internal matrix, exclude nx (end of matrix)    
if n<nx
    
beta_1 = (dx(n)^2/12)*(b_i(n)/k_i(n));
beta_2 = (dx(n+1)^2/12)*(b_i(n+1)/k_i(n+1));

K_11 = (k_i(n)/dx(n)*(1/(1+beta_1)) + (dx(n)/2)*b_i(n)); 

K_22 = (k_i(n+1)/dx(n+1))*(1/(1+beta_2)) + (dx(n+1)/2)*b_i(n+1);

K_21 = (-k_i(n)/dx(n))*(1/(1+beta_1)); 

K_12 = (-k_i(n+1)/dx(n+1))*(1/(1+beta_2)); 

K_i = [K_11,K_12;K_21,K_22]; 

    A(n,n) = K_i(1,1) + K_i(2,2); 
    A(n,n+1) = K_i(1,2);
    
%Prevents generation of i-1 in first row 
    if n>1
    A(n,n-1) = K_i(2,1);
    
    end
    
%Generate end values of matrix, excludes nx+1 value 
else
    
beta_2 = (dx(n)^2/12)*(b_i(n)/k_i(n));
    
K_21 = (-k_i(n)/dx(n))*(1/(1+beta_2)); 

K_22 = (k_i(n)/dx(n))*(1/(1+beta_2)) + (dx(n)/2)*b_i(n);

%Add to end of Global Matrix 
A(n,n) = K_22;
A(n,n-1) = K_21;
    
end

end

end