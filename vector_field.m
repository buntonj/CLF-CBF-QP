% GENERATE VECTOR FIELD FOR INPUTS
xrange = [-4:0.1:4];
yrange = [0:0.1:8];
[plotx, ploty] = meshgrid(xrange,yrange);
plotux = plotx;
plotuy = ploty;
[Q, q, r] = get_ellipse_parameters();

if method == 'standard'
    for i = 1:size(plotx,1)
        for j = 1:size(plotx,2)
            state = [plotx(i,j) ploty(i,j)];
            if state == q
                plotux(i,j) = 0;
                plotuy(i,j) = 0;
                plotx(i,j) = 0;
                ploty(i,j) = 0;
                continue
            end

            u_opt = CLF_CBF_QP(state);
            plotux(i,j) = u_opt(1)/sqrt(u_opt(1)^2 + u_opt(2)^2);
            plotuy(i,j) = u_opt(2)/sqrt(u_opt(1)^2 + u_opt(2)^2);
            %for not-normalized vector field
            %plotux(i,j) = u_opt(1);
            %plotuy(i,j) = u_opt(2);
        end
    end
elseif method == 'jankovic'
        for i = 1:size(plotx,1)
        for j = 1:size(plotx,2)
            state = [plotx(i,j) ploty(i,j)];
            if state == q
                plotux(i,j) = 0;
                plotuy(i,j) = 0;
                plotx(i,j) = 0;
                ploty(i,j) = 0;
                continue
            end

            u_opt = Jankovic_CLF_CBF_QP(state,Q,q,r);
            plotux(i,j) = u_opt(1)/sqrt(u_opt(1)^2 + u_opt(2)^2);
            plotuy(i,j) = u_opt(2)/sqrt(u_opt(1)^2 + u_opt(2)^2);
            %for not-normalized vector field
            %plotux(i,j) = u_opt(1);
            %plotuy(i,j) = u_opt(2);
        end
    end
end