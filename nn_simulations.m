%In this file, we consider a p by p array of cells with nearest and
%next-nearest neighbor couplings and periodic boundary conditions. We
%assume that the cells communicate via Notch signaling and simulate their
%cell-fates, modifying the model from Williamson et al., attached in the
%repository (williamson-model.pdf).
clear
 
%size of array of cells
p=16;


%%Construct adjacency matrix
%initialize adjacency matrix
A=zeros(p^2,p^2);

%coupling strengths of nearest and next-nearest neighbors
near=3;
n_near=1;

%Iterating through every cell in the array
%Loop from 0 to p-1 to match index with modular arithmetic
for i=0:p-1
    for j=0:p-1
    
    %See in_conv docstring
    %Converting (i,j) coordinate in cell array to cell # in ordering from
    %left to right and top to bottom
    row=in_conv(i+1,j+1,p);

    %inputs into the node of interest
    in1=in_conv(mod(i-1,p)+1,mod(j-1,p)+1,p);
    in2=in_conv(mod(i-1,p)+1,j+1,p);
    in3=in_conv(mod(i-1,p)+1,mod(j+1,p)+1,p);
    in4=in_conv(i+1,mod(j-1,p)+1,p);
    in5=in_conv(i+1,mod(j+1,p)+1,p);
    in6=in_conv(mod(i+1,p)+1,mod(j-1,p)+1,p);
    in7=in_conv(mod(i+1,p)+1,j+1,p);
    in8=in_conv(mod(i+1,p)+1,mod(j+1,p)+1,p);
    
    %nearest neighbors
    A(row,in2)=A(row,in2)+near;
    A(row,in4)=A(row,in4)+near;
    A(row,in5)=A(row,in5)+near;
    A(row,in7)=A(row,in7)+near;
    
    %next nearest neighbors
    A(row,in1)=A(row,in1)+n_near;
    A(row,in3)=A(row,in3)+n_near;
    A(row,in6)=A(row,in6)+n_near;
    A(row,in8)=A(row,in8)+n_near;
    end
end

%%Simulate ODE
%auxiliary functions f,g (Michaelis-Menten))
k=2;
h=2;
alpha=1; %Hill constants alpha, beta
beta=1;

f = @(x)  x^k/(alpha+x^k);
g = @(x)  1/(1+beta*x^h);

%other constants
gamma=1; %ratio of decay of M/N
mu=1.3; %ratio of EGF/Wnt mediated MapK activation

%Bifurcation parameters
E=7;
W=7;



%cells0 is a 2*p^2 x 1 matrix of initial conditions
%two vector entries represent one cell; one entry is the MapK (Delta)
%activity and the other is the Notch activity


%seed value
rng(100)
%Random initial conditions
cells0=rand(2*p^2,1);

%Time interval for ODE simulation
tspan = [0,500];

%Define vector field using vecfield function at the end of the file
dv = @(t,y) vecfield(y,E,W,f,g,mu,gamma,A);

%Solve ODE with built-in solver
[~,y] = ode45(dv,tspan,cells0);


%%Plot results
%Take final state of the system
yfinal = y(end,:);

%Separate yfinal vector into MapK and Notch components
MapK=zeros(length(yfinal)/2,1);
Notch=zeros(length(yfinal)/2,1);
for i=1:length(yfinal)/2
    Notch(i)=yfinal(2*i-1);
    MapK(i)=yfinal(2*i);

end

%reshape MapK and Notch vectors into an array for plotting in a heatmap
MapK=reshape(MapK,p,p)';
Notch=reshape(Notch,p,p)';

%Creating heatmaps
figure()
heatmap(Notch)
title("Notch")
colorbar("off")

figure()
heatmap(MapK)
title("MapK")
colorbar("off")


%%Defining the vector field for the array of cells


function Mbar=compute_Mbar(cells,A)
%Auxiliary function that averages the MapK from neighboring cells 
%Inputs:
%   cells: A vector of the Notch and MapK concentrations for every cell in
%   the array. It should follow the pattern (N1,M1,N2,M2,...)
%   A: The adjacency matrix of the array of cells
    
    %Compute the valence of each cell for averaging
    valence=A*ones(length(cells)/2,1);
    valence=valence(1);
    
    %Select MapK values from entire cell array (containing MapK and Notch)
    MapK=zeros(length(cells)/2,1);
    for i=1:length(cells)/2
        MapK(i)=cells(2*i);
    end
    
    %Multiplying the adjacency matrix A by the MapK vector sums all the
    %MapK in the neighbors. Average by multiplying by 1/valence.
    Mbar=(1/valence)*(A*MapK);
end
    

function dv = vecfield(cells,E,W,f,g,mu,gamma,A)
%Computes the vector field of the array of cells given the cells' current
%state and the auxiliary inputs below. For detailed information about the
%parameter meaning, see williamson-model.pdf and the supplementary
%information for their paper.
%Inputs
%   cells: A vector of the Notch and MapK concentrations for every cell in
%   the array. It should follow the pattern (N1,M1,N2,M2,...)
%   E,W: Scalar values that represent the bifurcation parameters
%   f,g: Auxiliary functions (previously defined)
%   mu,gamma: Additional scalar parameters
%   A: The adjacency matrix of the array of cells

    %Computing Mbar vector using previously defined function
    Mbar=compute_Mbar(cells,A);

    %Defining another parameter for the model
    lambda=(E+mu*W)/(mu+1);


    %Initializing vector field value
    dv = zeros(length(cells),1);
    %Loop through each cell
    for i=1:length(cells)/2
        
        %ODE Equations
        dv(2*i-1) = W*f(Mbar(i)) - cells(2*i-1);
        dv(2*i) = gamma*(lambda*g(cells(2*i-1))-cells(2*i));
    end
end


