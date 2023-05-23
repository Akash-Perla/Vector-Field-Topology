
%  Function to generate a sample 2D vector field
% doc randi

function [X,Y, ux, uy,ch] = test_vector_field()

Nt = 40;
xL = -3;
xU = 3;
yL = -2;
yU = 2;
x = linspace(xL,xU,Nt);
y = linspace(yL,yU,Nt);

[X,Y] = meshgrid(x,y);

dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);
[Ny,Nx] = size(X); 
[D2x,D2y,DD2x,DD2y] = F_diff_mat_2D_v2(Ny,Nx);  % Loading the finite difference matrices
ch = input("Choose a test field by selecting number between 1 and 4:");
fprintf("Trying Test Field %d\n",ch);

if( ch == 1)
    % Test field 1
    ux = Y.^2 - X.^2;
    uy = X.^2 + Y.^2 - 2;
elseif( ch == 2)
    % Test field 2
    ux = X.^2 - Y.^2;
    uy = X.^2 - sin(pi*Y.^2);
elseif( ch == 3)
    % Test field 3
    ux = Y.^2 - sin(0.5*pi*X.^2);
    uy = X.^2 + Y.^2 - 2;
else
    % Test field 4
    ux = Y.^2 - cos(0.5*pi*X.^2);
    uy= X.^2 + Y.^2 - 2;
end





