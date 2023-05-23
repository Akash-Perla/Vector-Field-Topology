
% Simplification of vector Field Topology (2D)


clear all;
close all;
clc; clf;


% Analysis parameters
% ======================================================================>>>
N_init_guess = 20;% Guess for number of critical points

min_tol = 2e-2;

rlim = 1000;     
                  
min_dist_tol = 2;

locus_dist_tol = 1e-2;   




%  Loading the vector field and defining the spatial grid


[X,Y,ux,uy,ch] = test_vector_field();     % Load a test vector field
fprintf("Displaying simplified vector field topology of test field %d\n",ch);

% Lower and upper bound of x and y
xL = min(X(:)); xU = max(X(:)); 
yL = min(Y(:)); yU = max(Y(:)); 

[Ny,Nx] = size(ux);     % Size of the vector field


u2 = ux.^2 + uy.^2;     % Squared norm function

quiver(X,Y,ux,uy);

dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);

% D2x = d/dx, DD2x = d^2/dx^2
[D2x,D2y,DD2x,DD2y] = F_diff_mat_2D_v2(Ny,Nx);  % Loading the finite difference matrices

J11 = reshape(D2x*ux(:)./(2*dx), size(X));  % dux/dx
J12 = reshape(D2y*ux(:)./(2*dy), size(X));  % dux/dy
J21 = reshape(D2x*uy(:)./(2*dx), size(X));  % duy/dx
J22 = reshape(D2y*uy(:)./(2*dy), size(X));  % dut/dy



fprintf('Critical point search:\n');
fprintf('==============================================\n');


counter = 0;

fprintf('Guesses of critical points (%i)\n', N_init_guess);
[u2s, ind_s] = sort(u2(:));

for m = 1:N_init_guess
    counter = counter + 1;
    xc_guess(counter) = X(ind_s(m));
    yc_guess(counter) = Y(ind_s(m));
%     fprintf('(%1.4f %1.4f) \n',xc_guess(counter),yc_guess(counter) );
end

u2min = min(u2(:));



fprintf('\nRevising\n');
fprintf('===================================================\n')
counter2=0;
for m = 1:length(xc_guess)
    temp = interp2(X,Y,u2,xc_guess(m),yc_guess(m));
    fprintf('(%1.4f,%1.4f) -> diff. from minimum = %1.4f', xc_guess(m),yc_guess(m),abs(temp-u2min));
    if abs(temp-u2min) <=  min_tol & xc_guess(m)^2 + yc_guess(m)^2 < rlim^2
        fprintf(': Accepted');
        counter2 = counter2 + 1;
        xc_guess_v(counter2) = xc_guess(m);
        yc_guess_v(counter2) = yc_guess(m);
    else
        fprintf(': Rejected');
    end
    fprintf('\n');
end

fprintf('\nRevised no. of critical points = %d\n' ,counter2);
Mtemp = [xc_guess_v' yc_guess_v'];
fprintf('(%1.4f %1.4f) \n',Mtemp' );
fprintf('===================================================\n')

counter = 0;
for n_critical = 1:length(xc_guess_v)
  
    x0 = xc_guess_v(n_critical);
    y0 = yc_guess_v(n_critical);
    
    fprintf('\nSimplification of critical point %i (%1.3f,%1.3f) \n' ,n_critical, x0, y0);

    
    xc_temp = x0; yc_temp = y0;
%     [xc_temp,yc_temp] = GridSearch(X,Y,ux,uy,x0,y0);
    [xc_temp,yc_temp,err_flag] = NewtonRaphson2D(X,Y,ux,uy,J11,J12,J21,J22,xc_temp,yc_temp);  % current critical point
    
    
    % Check if the current critical point is very close to a previously
    % calculated critical point. In such a case, set flag = 1 and do not
    % consider the current point as a new critical point.
    min_dist = inf;
    
    % calculate distance from the current critical point to all other
    % critical points. Then find the minimum distance.
    for q = 1:counter                                                      % Loop over all previous critical points
        dist = sqrt( (xc_temp - xc(q)).^2 + (yc_temp - yc(q)).^2 );        % Distance 
        if dist < min_dist
            min_dist = dist;
        end
    end
    
    if min_dist < min_dist_tol*dx
        err_flag = 1;
        fprintf('Duplicate found. Separation = %1.4f \n',min_dist);
    end
    
    if err_flag ~= 1
        
        counter = counter + 1;
        xc(counter) = xc_temp;
        yc(counter) = yc_temp;
    end
    
    fprintf('Refined root: (%1.4f, %1.4f)\n\n', xc(counter),yc(counter));
end
plotting_topology(X,Y,ux,uy,xc,yc);



        








