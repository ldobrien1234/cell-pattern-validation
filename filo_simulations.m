%In this file, we consider a p by p array of cells long-range couplings via
% filopodia and periodic boundary conditions. We assume that the cells 
% communicate via Notch signaling and simulate their cell-fates, modifying 
% the model from Williamson et al., attached in the repository 
%(williamson-model.pdf).
clear

%size of array of cells
p=50;

%%Construct adjacency matrix
%initialize adjacency matrix
A=zeros(p^2,p^2);

%Iterating through every cell in the array
%Loop through 0 to p-1 to match index with modular arithmetic
for i=0:p-1
    for j=0:p-1

    %See in_conv.m docstring 
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

    %filopodia
    fil1=in_conv(i+1,mod(j-2,p)+1,p);
    fil1u=in_conv(mod(i+1,p)+1,mod(j-2,p)+1,p);
    fil1d=in_conv(mod(i-1,p)+1,mod(j-2,p)+1,p);
    fil2=in_conv(i+1,mod(j-3,p)+1,p);
    fil2u=in_conv(mod(i+1,p)+1,mod(j-3,p)+1,p);
    fil2d=in_conv(mod(i-1,p)+1,mod(j-3,p)+1,p);
    fil3=in_conv(i+1,mod(j-4,p)+1,p);
    fil3u=in_conv(mod(i+1,p)+1,mod(j-4,p)+1,p);
    fil3d=in_conv(mod(i-1,p)+1,mod(j-4,p)+1,p);

    fil7=in_conv(i+1,mod(j+2,p)+1,p);
    fil7u=in_conv(mod(i+1,p)+1,mod(j+2,p)+1,p);
    fil7d=in_conv(mod(i-1,p)+1,mod(j+2,p)+1,p);
    fil8=in_conv(i+1,mod(j+3,p)+1,p);
    fil8u=in_conv(mod(i+1,p)+1,mod(j+3,p)+1,p);
    fil8d=in_conv(mod(i-1,p)+1,mod(j+3,p)+1,p);
    fil9=in_conv(i+1,mod(j+4,p)+1,p);
    fil9u=in_conv(mod(i+1,p)+1,mod(j+4,p)+1,p);
    fil9d=in_conv(mod(i-1,p)+1,mod(j+4,p)+1,p);



    %nearest neighbors
    A(row,in2)=A(row,in2)+1;
    A(row,in4)=A(row,in4)+1;
    A(row,in5)=A(row,in5)+1;
    A(row,in7)=A(row,in7)+1;
    
    %next nearest neighbors
    A(row,in1)=A(row,in1)+1;
    A(row,in3)=A(row,in3)+1;
    A(row,in6)=A(row,in6)+1;
    A(row,in8)=A(row,in8)+1;
    
    %filopodia
    A(row,fil1)=A(row,fil1)+3;
    A(row,fil1u)=A(row,fil1u)+3;
    A(row,fil1d)=A(row,fil1d)+3;
    A(row,fil2)=A(row,fil2)+3;
    A(row,fil2u)=A(row,fil2u)+3;
    A(row,fil2d)=A(row,fil2d)+3;
    A(row,fil3)=A(row,fil3)+3;
    A(row,fil3u)=A(row,fil3u)+3;
    A(row,fil3d)=A(row,fil3d)+3;


    A(row,fil7)=A(row,fil7)+3;
    A(row,fil7u)=A(row,fil7u)+3;
    A(row,fil7d)=A(row,fil7d)+3;
    A(row,fil8)=A(row,fil8)+3;
    A(row,fil8u)=A(row,fil8u)+3;
    A(row,fil8d)=A(row,fil8d)+3;
    A(row,fil9)=A(row,fil9)+3;
    A(row,fil9u)=A(row,fil9u)+3;
    A(row,fil9d)=A(row,fil9d)+3;
    end
end


%%Simulating ODE

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

%Bifurcation Parameters
E=6;
W=6;

%cells0 is a 2*p^2 x 1 matrix of initial conditions
%two vector entries represent one cell; one entry is the MapK
%activity and the other is the Notch activity

%seed value
rng(100)

%initial condition
cells0=randn(1,2*p^2)*10^(-3);

%Time interval for ODE simulations
tspan = [0,1200];

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
h1=heatmap(Notch);
min_val = min(Notch,[],"all");
max_val = max(Notch,[],"all");
h1.ColorLimits = [min_val(1) , max_val(1)]; 

% Define custom colormap using hex codes 
lowColor = hex2rgb('#FFFFFF');
highColor = hex2rgb('#FA8775');
rspan = linspace(lowColor(1),highColor(1))';
gspan = linspace(lowColor(2),highColor(2))';
bspan = linspace(lowColor(3),highColor(3))';


customColormap = colormap([rspan,gspan,bspan]);
% Apply custom colormap 
h1.Colormap = customColormap;
%h1.GridVisible = "off";

%Changing the color of the gridlines
h1_struct = struct(h1).Heatmap;
h1_grid = struct(h1_struct).Grid;
h1_grid.ColorData = uint8([200;200;200;125]);


figure()
h2=heatmap(MapK);
min_val = min(MapK,[],"all");
max_val = max(MapK,[],"all");
h2.ColorLimits = [min_val(1) , max_val(1)]; 

% Define custom colormap using hex codes 
lowColor = hex2rgb('#FFFFFF');
highColor = hex2rgb('#FA8775');
rspan = linspace(lowColor(1),highColor(1))';
gspan = linspace(lowColor(2),highColor(2))';
bspan = linspace(lowColor(3),highColor(3))';


customColormap = colormap([rspan,gspan,bspan]);
% Apply custom colormap 
h2.Colormap = customColormap; 
%h2.GridVisible = "off";

%Changing the color of the gridlines
h2_struct = struct(h2).Heatmap;
h2_grid = struct(h2_struct).Grid;
h2_grid.ColorData = uint8([200;200;200;125]);

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




