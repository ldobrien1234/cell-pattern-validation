%In this file, we consider a p by p array of cells with nearest and
%next-nearest neighbor couplings and periodic boundary conditions. We
%compute the eigenvector corresponding to the minimal eigenvalue of the
%adjacency matrix, which dictates the preferred pattern of the system.

clear


%dimension of cell array
p=16;

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

%Built-in function to compute the eigenvalues D and eigenvectors V
[V,D]=eig(A);

%Sort eigenvalues to ensure we consider the smallest one
d=sort(diag(D));


%Create heat array to visualize minimal eigenvalue
figure()
heat_array=reshape(V(:,1),[p,p])';
h=heatmap(heat_array);
% Set color limits
min_val = min(heat_array);
max_val = max(heat_array);
h.ColorLimits = [min_val(1) , max_val(1)]; 

% Define custom colormap using hex codes 
lowColor = hex2rgb('#FFFFFF');
highColor = hex2rgb('#FA8775');
rspan = linspace(lowColor(1),highColor(1))';
gspan = linspace(lowColor(2),highColor(2))';
bspan = linspace(lowColor(3),highColor(3))';


customColormap = colormap([rspan,gspan,bspan]);
% Apply custom colormap 
h.Colormap = customColormap; 

%Adjust the color of the gridlines
h_struct = struct(h).Heatmap;
h_grid = struct(h_struct).Grid;
h_grid.ColorData = uint8([200;200;200;125]);
colorbar("off")
