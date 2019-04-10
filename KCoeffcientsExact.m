function [A] = KCoeffcientsExact(dx,b_i,k_i,nx)

alpha_1 = sqrt(b_i(1)/k_i(1)); 


K_11 = k_i(1)*(alpha_1/(sinh(alpha_1*dx(1)))) * (cosh(alpha_1*dx(1)));
K_12 = k_i(1)*(alpha_1/(sinh(alpha_1*dx(1))));



A(1,1) = K_11;
A(1,2) = -K_12;

for i =2:nx-1
    
    alpha_1 = sqrt(b_i(i)/k_i(i));
    alpha_2 = sqrt(b_i(i+1)/k_i(i+1));


K_11 = k_i(i)*(alpha_1/(sinh(alpha_1*dx(i)))) * (cosh(alpha_1*dx(i)));
K_12 = -k_i(i)*(alpha_1/(sinh(alpha_1*dx(i))));

K_22 = k_i(i+1)*(alpha_2/(sinh(alpha_2*dx(i+1)))) * (cosh(alpha_2*dx(i+1)));
K_21 = -k_i(i+1)*(alpha_2/(sinh(alpha_2*dx(i+1))));

A(i,i-1) = K_12;

A(i,i) = K_11+K_22;
A(i,i+1) = K_21;


end

alpha_2 = sqrt(b_i(nx)/k_i(nx));

K_22 = k_i(nx)*(alpha_2/(sinh(alpha_2*dx(nx)))) * (cosh(alpha_2*dx(nx)));
K_21 = -k_i(nx)*(alpha_2/(sinh(alpha_2*dx(nx))));

A(nx,nx-1) = K_21;
A(nx,nx) = K_22;

