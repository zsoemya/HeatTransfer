function [A] = GlobalMatrix4thFDM(dx,b_i,k_i,nx)
%Generates Global Matrix A without Boundry Conditions for 4th Order FDM
%   Detailed explanation goes here


%Generate Initial Line of Coeffcients for A
k_intial = [1,k_i(1)];
b_intial = [0,b_i(1)];

K_initial = KCoeffcients4thFDM(dx,b_intial,k_intial,1); 

A(1,1) = K_initial(1,1) + K_initial(2,2);
A(1,2) = K_initial(1,2); 

%Add Additional Length to b and k coeffcient matrix to generate last term
%U(n+1)
b_i(nx+1,1) = 0;
k_i(nx+1,1) = 1;
dx(1,nx+1) = dx(1,nx); 


for i=2:nx
    
    K_i = KCoeffcients4thFDM(dx,b_i,k_i,i);
    
    A(i,i) = K_i(1,1) + K_i(2,2); 
    A(i,i-1) = K_i(2,1);
    
    %Prevent generation of extra column at end of domain 
    if i<nx
    A(i,i+1) = K_i(1,2);
    end   
end

end



