clear
clc

%Problem Initialisation 
%Allows user to input all required varibles to solve problem.
%To avoid excessive coding minimal error returns have been implemented and
%common sense should be used, ie. do not input a negative delta x value 

%%
%Domain User Inputs 

Start_X = 'x domain min value? ';
x_min = input(Start_X); 

Refine = 'Number of Mesh Refinements? ';
Iter = input(Refine); 


%Element User Inputs
Elements = 'Number of Material Elements? ';

No_Elements = input(Elements);


%Coeffcient and Element Length User Inputs 

for i = 1:No_Elements 
    
    x_El = 0; 
    
    K_Length = ['Length of Element ',num2str(i),'? (must be mutiple of Delta x) '];
    k_l(i) = input(K_Length);
    
    if i==1
    Total_Length(i) = k_l(i) + x_min;
    else
    Total_Length(i) = k_l(i) + Total_Length(i-1);
    end
    
    El_Delta = ['Delta x Value for Element: '];
    el_dx(i) = input(El_Delta);
    
    K_Value = ['k-Coeffcient of Element ',num2str(i),'? '];
    k(i) = input(K_Value);
    
    B_Value = ['b-Coeffcient of Element ',num2str(i),'? '];
    b(i) = input(B_Value);  

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Boundry Conditions User Inputs
%End Condtion 
End_BC_prompt = 'Specify End Boundry Condition [Y/N] ? ';

End_BC = input(End_BC_prompt,'s');

End_BC_Type = 0;

if End_BC == 'Y'
    
End_BC_Type_Prompt = 'BC Type: Flux BC =1, Essential BC = 2, Newtown BC = 3, Input Selection: ';

End_BC_Type = input(End_BC_Type_Prompt);

if End_BC_Type == 1
    
    Final_Q_Prompt = 'Specifiy Q(L) Value:';
    B_Q_L = input(Final_Q_Prompt); 
    
elseif End_BC_Type == 2
    
    Final_U_Prompt = 'Specify U(L) Value: ';
    
    T_L = input(Final_U_Prompt); 
    
    b_L = 10^20; 
    
elseif End_BC_Type == 3
    
    Beta_Prompt = 'Specify b(L) Value: ';
    
    b_L = input(Beta_Prompt); 
    
    U_L_Prompt = 'Specify T(L) Value: ';
    
    T_L = input(U_L_Prompt); 
    
end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Start Condtion 
Start_BC_prompt = 'Specify Start Boundry Condition [Y/N] ? ';

Start_BC_Type = 0;

Start_BC = input(Start_BC_prompt,'s');

if Start_BC == 'Y'
    
Start_BC_Type_Prompt = 'BC Type: Flux BC =1, Essential BC = 2, Newtown BC = 3, Input Selection: ';

Start_BC_Type = input(Start_BC_Type_Prompt);

if Start_BC_Type == 1
    
    Start_Q_Prompt = 'Specifiy Q(0) Value: ';
    B_Q_O = input(Start_Q_Prompt); 
    
elseif Start_BC_Type == 2
    

    Start_U_Prompt = 'Specify T(0) Value: ';
    
    T_0 = input(Start_U_Prompt); 
    
    b_0 = 10^20;
    
elseif Start_BC_Type == 3
    
    
    Beta_Prompt = 'Specify b(0) Value: ';
    
    b_0 = input(Beta_Prompt); 
    
    U_0_Prompt = 'Specify T(0) Value: ';
    
    T_0 = input(U_0_Prompt);
    
    
else
    
    return
end
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m = 1:Iter
    
    if m > 1
       el_dx = el_dx./2;
       x = 0;
    end

%Domain Construction 
    for i = 1:No_Elements

    if i==1
    x_El = x_min:el_dx(1):Total_Length(1);
    nx_El(1) = length(x_El);
    x(1,i:nx_El) = x_El;
    dx(1,i:nx_El) = el_dx(i);
    else
    x_El = Total_Length(i-1)+el_dx(i):el_dx(i):Total_Length(i);
    nx_El(i) = length(x_El);
    nx = length(x);
    x(1,(nx+1):(nx_El(i)+nx)) = x_El;
    dx(1,(nx+1):(nx_El(i)+nx)) = el_dx(i);
    end
    
    if i ==1
    length_index(i) = (k_l(i)/el_dx(i)) + 1;
    else
    length_index(i) = (k_l(i)/el_dx(i)) + length_index(i-1); 
    end
    
    if i == 1
    
    k_i(1:length_index(i),1) = k(i);
    b_i(1:length_index(i),1) = b(i);
    
    else
       
    k_i(length_index(i-1)+1:length_index(i),1) = k(i);
    b_i(length_index(i-1)+1:length_index(i),1) = b(i);  
        
    end
    
    end
    
%Add End Condtions     
nx = length(x); 

%Intialise B Matrix 
B = zeros(nx,1);

%End Condtions 
if End_BC_Type == 1
    
    B(nx,1) = B_Q_L;
    
elseif End_BC_Type == 2
    
    B(nx,1) = T_L*b_L;
    
elseif End_BC_Type == 3
    
    B(nx,1) = T_L*b_L; 
end

%Start Conditions 
if Start_BC_Type == 1
    
    B(1,1) = B_Q_O;
    
elseif Start_BC_Type == 2
    
    B(1,1) = T_0*b_0; 
    
elseif Start_BC_Type == 3

    B(1,1) = T_0*b_0; 
    
end

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%



%%

%Mode Settings for Global Matrix 
FDM_2nd = 1;
FEM_p1 = 2;

%Global Matrix Generation 
A_FDM2 = GlobalMatrix(dx,b_i,k_i,nx,FDM_2nd);

A_FEM = GlobalMatrix(dx,b_i,k_i,nx,FEM_p1);

A_FDM4 = KCoeffcients4thFDM(dx,b_i,k_i,nx);

A_Exact = KCoeffcientsExact(dx,b_i,k_i,nx);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%Boundry Condtion Modification

%End Condtion Modification
if End_BC_Type == 2 || End_BC_Type == 3
%Newtownian End Conditions
A_FDM2(nx,nx) = A_FDM2(nx,nx) + b_L;
A_FEM(nx,nx) = A_FEM(nx,nx) + b_L;
A_FDM4(nx,nx) = A_FDM4(nx,nx) + b_L;
A_Exact(nx,nx) = A_Exact(nx,nx) + b_L;

end

if Start_BC_Type == 2 || Start_BC_Type == 3
    
%Newtownian End Conditions
A_FDM2(1,1) = A_FDM2(1,1) + b_0;
A_FEM(1,1) = A_FEM(1,1) + b_0;
A_FDM4(1,1) = A_FDM4(1,1) + b_0;
A_Exact(1,1) = A_Exact(1,1) + b_0;
    
end



%%
%Solving 
U_FDM2 = linsolve(A_FDM2,B);
U_FEM = linsolve(A_FEM,B);
U_FDM4 = linsolve(A_FDM4,B);
U_Exact = linsolve(A_Exact,B);

disp('Ignore Matrix Scaling Error Above')

%Caluclate Q (Heat Flux)
Q_FDM2 = Q_Flux(k_i,b_i,dx,U_FDM2);
Q_FEM = Q_Flux(k_i,b_i,dx,U_FEM);
Q_FDM4 = Q_Flux(k_i,b_i,dx,U_FDM4);
Q_Exact = Q_Flux(k_i,b_i,dx,U_Exact);


%Store Results for Convergence  
dx_Plot(m,1) = max(dx);

FDM2_Center(m,1) = U_FDM2(ceil(nx/2),1);
FEM_Center(m,1) = U_FEM(ceil(nx/2),1);
FDM4_Center(m,1) = U_FDM4(ceil(nx/2),1);
Exact_Center(m,1) = U_Exact(ceil(nx/2),1);

end

%%
%Results

FDM2_Error = relError(dx_Plot,FDM2_Center,Exact_Center);
FEM_Error = relError(dx_Plot,FEM_Center,Exact_Center);
FDM4_Error = relError(dx_Plot,FDM4_Center,Exact_Center);

%Plots 

%Plot of highest refined mesh for each method 
figure(1)
hold on
plot(x,U_Exact,'Linewidth',1)
plot(x,U_FDM2)
plot(x,U_FEM)
plot(x,U_FDM4)
xlabel('Location x (cm)');
ylabel('Temperture T(x) (C)');
grid on


legend({'Exact','2nd Order FDM','p=1 FEM','4th Order FDM'},'Location','northwest');

%Table Data 
X = transpose(x); 
Exact = U_Exact;
SecondOrderFDM = U_FDM2;
p1FEM = U_FEM;
FourthOrderFDM = U_FDM4;

Q_SecondOrderFDM = Q_FDM2;
Q_p1FEM = Q_FEM;
Q_FourthOrderFDM = Q_FDM4;

%Errors
SecondOrderFDMError = abs(U_Exact-SecondOrderFDM);
p1FEMError = abs(U_Exact-p1FEM);
FourthOrderFDMError = abs(U_Exact-FourthOrderFDM);


%Produce Table of Results for highest refined mesh U results 
T = table(X,U_Exact,SecondOrderFDM,SecondOrderFDMError,p1FEM,p1FEMError,FourthOrderFDM,FourthOrderFDMError);

disp(T);

%Produce Table of Results for highest refined mesh Q results 
T_Q = table(X,Q_Exact,Q_SecondOrderFDM,Q_p1FEM,Q_FourthOrderFDM);

disp(T_Q);

%Plot convergence graph
figure(2)
grid on
relerror = loglog(dx_Plot,FDM2_Error,'-*','Linewidth',2);
hold on
loglog(dx_Plot,FEM_Error,'-*','Linewidth',2);
loglog(dx_Plot,FDM4_Error,'-*','Linewidth',2);
set(gca,'Xdir','reverse')
xlabel('Log|\Deltax|');
ylabel('Log|Relative Error|(%)');
legend('2nd Order FDM','p=1 FEM','4th Order FDM')
grid on

%%


