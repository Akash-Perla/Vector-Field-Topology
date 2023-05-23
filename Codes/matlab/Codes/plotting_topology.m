function []=plotting_topology(X,Y,ux,uy,xc,yc)

xL = min(X(:)); xU = max(X(:)); 
yL = min(Y(:)); yU = max(Y(:)); 


dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);
[Ny,Nx] = size(ux); 
% D2x = d/dx, DD2x = d^2/dx^2
[D2x,D2y,DD2x,DD2y] = F_diff_mat_2D_v2(Ny,Nx);

J11 = reshape(D2x*ux(:)./(2*dx), size(X));  % dux/dx
J12 = reshape(D2y*ux(:)./(2*dy), size(X));  % dux/dy
J21 = reshape(D2x*uy(:)./(2*dx), size(X));  % duy/dx
J22 = reshape(D2y*uy(:)./(2*dy), size(X));  % dut/dy

fprintf('Final no. of critical points found = %d \n',length(xc));
Mtemp = [xc' yc'];
% fprintf('(%1.4f %1.4f) \n',Mtemp' );

clr_arr = [0.0   0.4   0.8;
    1.0   0.0   0.0;
    0.0   0.6   0.0;
    1.0   1.0   0.0;
    0.5   0.5   0.5;
    0.8   0.0   0.8;
    ];
cr_clss = zeros(1,size(xc,1));
txt_arr = {'Att. focus', 'Rep. focus', 'Saddle' 'Center' ,'Att. node','Rep. node'};


for m = 1:length(xc)
    
    temp_J = JacobianInterp(J11,J12,J21,J22,X,Y,xc(m),yc(m));  % Calculating the Jacobian matrix at critical point
    [temp_V,temp_D] = eig(temp_J);  % Calculating the Eigen values and Eigen vectors
    
    temp_lambda = diag(temp_D);     % Eigen values
    
    cr_clss(m) = critical_class(temp_lambda);  % Classifying the critical points
    
    figure(1); hold on;
    plot(xc(m),yc(m),'Marker','o','MarkerFaceColor',clr_arr(cr_clss(m),:),'markersize',8,'color','none');
    axis([xL xU yL yU]);
    
    daspect([1 1 1]);
    set(gca,'Position',[0.1 0.1 .65 .65]);
end


xsink = xc(cr_clss == 1 | cr_clss == 5);
ysink = yc(cr_clss == 1 | cr_clss == 5);

cnt = 0;
for m = 1:length(xc)
    
    % Check critical point class. For sinks, the direction has to be
    % reversed (the lines propagate in the opposite direction of the field
    % vector).
    if cr_clss(m) == 1 | cr_clss(m) == 4 | cr_clss(m) == 5
        drct = -1;
    else
        drct = 1;
    end
    
    
    x0 = xc(m);
    y0 = yc(m);
    J0 = JacobianInterp(J11,J12,J21,J22,X,Y,x0,y0);
    [temp_V,temp_D] = eig(temp_J);
    
    
    temp_V = real(temp_V);    
    
    
    
    if rank(temp_V) == 2
        
        xslope_set = [temp_V(1,1) -temp_V(1,1) temp_V(1,2) -temp_V(1,2)];
        yslope_set = [temp_V(2,1) -temp_V(2,1) temp_V(2,2) -temp_V(2,2)];
    else
        xslope_set = [temp_V(1,1) -temp_V(1,1) temp_V(2,1) -temp_V(2,1)];
        yslope_set = [temp_V(2,1) -temp_V(2,1) -temp_V(1,1) temp_V(1,1)];
    end
    
    
    
    for eg_vct = 1:length(xslope_set)
        xslope = xslope_set(eg_vct);
        yslope = yslope_set(eg_vct);
        
        
        
        xy_locus = int_locus(X,Y,ux,uy,x0,y0,xslope,yslope,xsink,ysink,drct);
        
        cnt = cnt + 1;
        xlocus{cnt} = xy_locus(:,1);
        ylocus{cnt} = xy_locus(:,2);
        cr_clss_save(cnt) = cr_clss(m);
        
        % For saddle points, draw in +drct and -drct directions. There are
        % redundancies though. For each direction, there are just 2 lines,
        % instead of 4.
        if cr_clss(m) == 3
            xy_locus = int_locus(X,Y,ux,uy,x0,y0,xslope,yslope,xsink,ysink,-drct);
            cnt = cnt + 1;
            xlocus{cnt} = xy_locus(:,1);
            ylocus{cnt} = xy_locus(:,2);
            cr_clss_save(cnt) = cr_clss(m);
        end
        
    end
    
    fprintf('Finished calculating integral lines for critical point %d \n',m);
    
end




figure(1);
hold on;
for m = 1:length(xlocus)
    plot(xlocus{m},ylocus{m},'k');
    
    drawArrow(xlocus{m},ylocus{m},X,Y,ux,uy,clr_arr(cr_clss_save(m),:));
    
    
end

xlabel('x'); ylabel('y');
axis([xL xU yL yU]);

daspect([1 1 1]);
set(gca,'Position',[0.1 0.1 .65 .65]);


for m = 1:6
    annotation('ellipse',[.82 m/10+.05 .015 .02],'Facecolor',clr_arr(m,:),'color','none');
    annotation('textbox',[.83 m/10+.085 0 0],'string',txt_arr{m});
    
end
pbaspect([1 1 1]);
saveas(gcf, 'Simplified_Topology.png');

