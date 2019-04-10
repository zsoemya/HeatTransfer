function [K_i] = KCoeffcients2ndFDM(dx,b_i,k_i,n)

K_11 = (k_i(n)/dx(n)) + (dx(n)/2)*b_i(n); 

K_22 = (k_i(n+1)/dx(n+1)) + (dx(n+1)/2)*b_i(n+1);

K_21 = (-k_i(n)/dx(n)); 

K_12 = (-k_i(n+1)/dx(n+1)); 

K_i = [K_11,K_12;K_21,K_22]; 

end