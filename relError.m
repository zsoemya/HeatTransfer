function [relativeError] = relError(dx ,computed ,exact)
%Find relative error between a computed and an exact value

nx = length(dx); 

for i = 1:nx

relativeError(i,1) = abs(exact(i,1) - computed(i,1))./(exact(i,1)*100);

end

end

