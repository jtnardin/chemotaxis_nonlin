%chemotaxis_computation written 11-7-17 by JTN to compute solution
%to chemotaxis equation

function u = chemotaxis_computation(q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A_pos,...
    A_neg,D_nonlin,D_lin,D_nonlin_bd,D_lin_bd)


%     q = [D,chi,n,alpha];

    v_ind_int = x_int + xn;
    v_xbd_0 = xbd_0 + xn;
    v_xbd_1 = xbd_1 + xn;
    
    %evaluate V
    Vx_c = @(v) q(2)*[diff(v);v(end)-v(end-1)]./(dx*(1+q(4)*v.^2));

    %crank nicholson
    theta = 0.5;

   
    %sigma for flux limiters
    sigma = @(r) (r+abs(r))./(1+abs(r));


    %initialize u
    u = zeros(2*xn,tn);
    u(:,1) = IC;
   
    sparse_id = sparse([xn+1:2*xn xn+1:2*xn],[1:xn xn+1:2*xn],...
        [ones(1,xn) -ones(1,xn)],2*xn,2*xn);


    for i = 2:tn

        
        v = dt/dx*Vx_c(u(xn+1:end,i-1));
        
        %get Ax matrices
        [Axpos,Axneg] = AMatrixCompute(u(:,i-1),v,x_int,xn,A_pos,A_neg);

        u(:,i) = (speye(2*xn) + theta*(Axpos+Axneg))\(speye(2*xn) - ...
            (1-theta)*(Axpos+Axneg) + D_nonlin(u,x_int) +...
            D_nonlin_bd(u,xbd_0,xbd_1) + D_lin(v_ind_int) + ...
            D_lin_bd(v_xbd_0,v_xbd_1) + sparse_id)*u(:,i-1);

        
    end

    
end
