%MLE_cost_D0 written 6-13-17 by JTN to do 29-point WLS scheme

function [J,OLS_SV,res,model] = MLE_cost(cell_data,...
    q_est,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A_pos,...
        A_neg,D_nonlin,D_lin,D_nonlin_bd,D_lin_bd,x,t,xdata,tdata)

    %update q_est vector based on sim_type
    
    alpha = 0;
    q_est = [q_est(1:2) ; alpha ; q_est(3:end)];
    
   %run simulation
    [model] = chemotaxis_computation(q_est,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A_pos,...
        A_neg,D_nonlin,D_lin,D_nonlin_bd,D_lin_bd,x,t,xdata,tdata);
    
    
    %each row of model corresponds to solution at given time point. In this
    %stat model, we only care about rows when t is a multiple of toCompare
    %same with data;
    
%     model = model(mod(tdata,toCompare)==0,:);
%     cell_data = cell_data(mod(tdata,toCompare)==0,:);
    
    %total data points considered
    N = numel(cell_data);
    
    
    %calculate residuals
    res = model(:) - cell_data(:);
    
    OLS_SV = 1/N*sum(res.^2); %slightly biased sample variance
    
    
    J = sum(res.^2);


end