%nonlin_taxis_fitting.m written  11-17-17 by JTN to perform
%OLS optimization for a nonlinear diffusion + chemotaxis equation in a 
%parallelized for loop to mean experimental data


function nonlin_taxis_fitting(m,pred_ind,simnum)

    %for bookkeeping
    well = m;
    
    %interp IC
    IC_type = 'interp';
        
            
    %load in data
    load('ind_cell_prof_data.mat')
    
    %helps with FRET interpolation
    dx_large = 10;
    
    tstart = 16;
    %choose data
    cell_data = avg_cell_data(ind_cell_data{m-1,2}(:,:,tstart:end),pred_ind)';
    FRET_data = avg_cell_data(ind_fret_data{m-1,2}(:,tstart:dx_large:end,tstart:end),pred_ind)';
    
        
    %initialize data grids
    [tndata,xndata] = size(cell_data);
    xdata = linspace(0,1,xndata);
    xdata2 = xdata(1:dx_large:end);
    tdata = 5:1/3:1/3*(tndata-1+5*3);
    %make ndgrid for interpolation
    [Xdata,Tdata] = ndgrid(xdata,tdata);
    [Xdata2,Tdata2] = ndgrid(xdata2,tdata);

    %generate grids for computation, also helpful for interpolation
    xn = 100;
    dt = 1e-3;
    [x,t] = grid_generate(xn,xdata(1),xdata(end),dt,tdata(1),tdata(end));
    tn = length(t);
    dx = x(2)-x(1);
    %interior, boundary points
    [x_int,xbd_0,xbd_1] = int_bd_def(xn);
    
    
    %%%% Now fit to migration data. First initialize q and cost vectors
    q_all = cell(simnum,1);
    q0_all = cell(simnum,1);
    J_all = zeros(simnum,1);

   

    %initial condition
    switch IC_type
    case 'interp'

            ICu = interp1(xdata,cell_data(1,:),x);
            cutoff_x = leading_edge_calc(ICu,x,.05,0);
            ICu(x>cutoff_x)=0;
            
            ICv = 1*ones(size(x));
            
            IC = [ICu(:) ; ICv(:)];

    case 'step'

            LE_loc = leading_edge_calc(smooth(cell_data(1,:)),xdata,0.8,0);
            IC = double(x<=LE_loc);
    end


    %boundary conditions
    BC_x_0 = @(t) 1;
    BC_x_1 = @(t) 0;


    [A_pos,~,~,A_neg,~,~,D_nonlin,D_nonlin_bd,D_lin,...
        D_lin_bd] = aMatrixConstruction(xn,dx);

    options = optimset('maxiter',100,'display','iter');
       

    %each row of data matrix corresponds to data at given time point. In this
    %stat model, we only care about rows when t is a multiple of toCompare

    %%% also only looking at every 5th point in space now (10-25)
    
    toCompare = 6;
    
    cell_data = cell_data(1:3*toCompare:end,1:5:end);

    tdata = tdata(1:3*toCompare:end);
    xdata = xdata(1:5:end);

    name_save = [num2str(toCompare)];
        
    
    for k = 1:simnum

        
         %q = [D, chi, k_v, D_v];
         q0_all{k} = [.0005 0.1 0.005 .01]'*2.*rand(4,1);
         
         LB = [zeros(4,1)];
         UB = [inf(4,1)];
         

        tic

        [q_all{k},J_all(k)] = fmincon(@(q) MLE_cost(cell_data,...
            q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A_pos,...
            A_neg,D_nonlin,D_lin,D_nonlin_bd,D_lin_bd,x,t,xdata,tdata),q0_all{k}...
            ,[],[],[],[],LB,UB,[],options);

        toc

                       
    end

    %save
    save(['/scratch/summit/jona8898/chem_fitting/taxis_fitting_well_' ...
        num2str(well) '_' name_save '_pred_'...
        num2str(pred_ind) '.mat' ],'q_all','q0_all','J_all')

end
