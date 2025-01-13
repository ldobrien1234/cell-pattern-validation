%In this file, we consider a p by p array of cells with ranged couplings 
% via filopodia and periodic boundary conditions. We compute the 
% eigenvector corresponding to the minimal eigenvalue of the adjacency 
% matrix, which dictates the preferred pattern of the system.


clear


%dimension of cell array
p=50;

%initialize adjacency matrix
A=zeros(p^2,p^2);

%Iterating through every cell in the array
%Loop from 0 to p-2 to match index with modular arithmetic
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

%Built-in function to compute the eigenvalues D and eigenvectors V
[V,D]=eig(A);

%Sort eigenvalues to ensure we consider the smallest one
d=sort(diag(D));

%Reshape two eigenvectors corresponding to minimum eigenvalue into an array
%for plotting
heat_array1=(reshape(V(:,1), [p,p]))';
heat_array2=(reshape(V(:,2),[p,p]))';

%Heatmap of first eigenvector
figure()
h1 = heatmap(heat_array1);
min_val = min(heat_array1,[],"all");
max_val = max(heat_array1,[],"all");
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

%Adjust the color of the gridlines
h1_struct = struct(h1).Heatmap;
h1_grid = struct(h1_struct).Grid;
h1_grid.ColorData = uint8([200;200;200;125]);

colorbar("off")


%Heatmap for second eigenvectors
figure()
h2 = heatmap(heat_array2);
min_val = min(heat_array2,[],"all");
max_val = max(heat_array2,[],"all");
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

%Adjust the color of the gridlines
h2_struct = struct(h2).Heatmap;
h2_grid = struct(h2_struct).Grid;
h2_grid.ColorData = uint8([200;200;200;125]);
colorbar("off")
