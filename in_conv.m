function n=in_conv(i,j,p)
%We have a p by p array of cells. To construct the adjacency matrix, we
%need an ordering of the cells 1,...,p^2. We choose to order them from left
%to right and top to bottom. This function converts the array index (i,j)
%to a number in the ordering given the dimension of the square array, p.
    n=p*(i-1)+j;
end
