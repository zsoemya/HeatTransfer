function [A] = GlobalMatrix(dx,b_i,k_i,nx,mode)
%Generates Global Matrix A without Boundry Conditions
%   Mode 1 is 2nd Order FDM
%   Mode 2 is p=1 FEM


%Generate Initial Line of Coeffcients for A
k_intial = [0,k_i(1)];
b_intial = [0,b_i(1)];

if mode == 1

K_initial = KCoeffcients2ndFDM(dx,b_intial,k_intial,1); 

elseif mode == 2
    
K_initial = KCoeffcientsp1FEM(dx,b_intial,k_intial,1);
    
end

A(1,1) = K_initial(1,1) + K_initial(2,2);
A(1,2) = K_initial(1,2); 

%Add Additional Length to b and k coeffcient matrix to generate last term
%U(n+1)
b_i(nx+1,1) = 0;
k_i(nx+1,1) = 0;
dx(1,nx+1) = dx(1,nx);


for i=2:nx
    
    if mode==1
    
    K_i = KCoeffcients2ndFDM(dx,b_i,k_i,i);
    
    elseif mode==2
     
    K_i = KCoeffcientsp1FEM(dx,b_i,k_i,i);  
        
    end
    
    A(i,i) = K_i(1,1) + K_i(2,2); 
    A(i,i-1) = K_i(2,1);
    
    %Prevent generation of extra column at end of domain 
    if i<nx
    A(i,i+1) = K_i(1,2);
    end   
end

end



