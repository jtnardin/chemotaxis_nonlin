%chemotaxis_computation written 11-7-17 by JTN to compute solution
%to chemotaxis equation

function udata = chemotaxis_computation(q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A_pos,...
    A_neg,D_nonlin,D_lin,D_nonlin_bd,D_lin_bd,x,t,xdata,tdata)


%     q = [D,chi,n,alpha];

    v_ind_int = x_int + xn;
    v_xbd_0 = xbd_0 + xn;
    v_xbd_1 = xbd_1 + xn;
    
    %evaluate V
    Vx_c = @(v) q(2)*[diff(v);v(end)-v(end-1)]./(dx*(1+q(3)*v.^2));

    %crank nicholson
    theta = 0.5;


    %initialize u
    u = zeros(2*xn,tn);
    u(:,1) = IC;
   
    

    for i = 2:tn

        
        sparse_id = sparse([xn+1:2*xn xn+1:2*xn],[1:xn xn+1:2*xn],...
         [zeros(1,xn) -q(4)*u(1:xn,i-1)'],2*xn,2*xn);   %,[q(5)*ones(1,xn) -ones(1,xn)]

        
        v = dt/dx*Vx_c(u(xn+1:end,i-1));
        
        %get Ax matrices
        [Axpos,Axneg] = AMatrixCompute(u(:,i-1),v,x_int,xn,A_pos,A_neg);

        u(:,i) = (speye(2*xn) + theta*(Axpos+Axneg))\(speye(2*xn) - ...
            (1-theta)*(Axpos+Axneg) + dt*q(1)*D_nonlin(u(:,i-1),x_int) +...
            dt*q(1)*D_nonlin_bd(u(:,i-1),xbd_1) + dt*q(5)*D_lin(v_ind_int) + ...
            dt*q(5)*D_lin_bd(v_xbd_0,v_xbd_1) + dt*sparse_id)*u(:,i-1);

        
    end
    
    
    [X,T] = meshgrid(x,t);
    
    udata = interp2(X,T,u(1:xn,:)',xdata,tdata');

    
end
