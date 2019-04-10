function Q = Q_Flux(k_i,b_i,dx,U)
%Caluculates the Q (Heat Flux) from U (Tempertures)

nx = length(U);

for v = 2:nx
    
Q(v,1) = -k_i(v)/dx(v)*U(v-1,1) + (k_i(v)/dx(v) + b_i(v)*(dx(v)/2))*U(v,1);

end


end

