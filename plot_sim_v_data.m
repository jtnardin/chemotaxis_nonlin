% plot_sim_v_data written 10-4-17 by JTN to simply take a parameter, q, and
% run the model against some data.

clear all; clc


%define parameters for simulation:
well = 3;
pred_ind = 1;

stat_model = 'fewer';
IC_type = 'interp';
fitpred = 'fit';

plot_res = 0;


%load parameter vector
% load(['taxis_fitting_well_' num2str(well) '_6_pred_' num2str(pred_ind) '.mat'])
q = [.0005 0.1 0.005 .01]'; %.05


alpha = 0;
q = [q(1:2) ; alpha ; q(3:end)];

%make sure q is a column vector
q = reshape(q,length(q),1);


%load in data
load('ind_cell_prof_data.mat')

%helps with FRET interpolation
dx_large = 10;
%choose data

if strcmp(fitpred,'fit')

    cell_data = avg_cell_data(ind_cell_data{well-1,2}(:,:,16:end),pred_ind)';
    FRET_data = avg_cell_data(ind_fret_data{well-1,2}(:,1:dx_large:end,16:end),pred_ind)';
    cell_data_std = ind_cell_data_sv{well-1,2}(:,16:end);

elseif strcmp(fitpred,'pred')


    cell_data = squeeze(ind_cell_data{well-1,2}(pred_ind,:,16:end))';
    FRET_data = squeeze(ind_fret_data{well-1,2}(pred_ind,1:dx_large:end,16:end))';
    cell_data_std = ind_cell_data_sv{well-1,2}(:,16:end);

else

    error('incorrect fitpred')

end

%initialize data grids
[tndata,xndata] = size(cell_data);
xdata = linspace(0,1,xndata);
xdata2 = xdata(1:dx_large:end);
tdata = 5:1/3:1/3*(tndata-1+5*3);
%make ndgrid for interpolation
[Xdata,Tdata] = ndgrid(xdata,tdata);
[Xdata2,Tdata2] = ndgrid(xdata2,tdata);


xn = 100;
dt = 1e-3;
[x,t] = grid_generate(xn,xdata(1),xdata(end),dt,tdata(1),tdata(end));
tn = length(t);
dx = x(2)-x(1);
%interior, boundary points
[x_int,xbd_0,xbd_1] = int_bd_def(xn);


%initial condition
switch IC_type
case 'interp'

        ICu = interp1(xdata,cell_data(1,:),x);
        cutoff_x = leading_edge_calc(ICu,x,.05,0);
        ICu(x>cutoff_x)=0;

%         ICu = double(x>.2).*double(x<=.5);

        ICv = 1*ones(size(x));

        IC = [ICu(:) ; ICv(:)];

    
case 'step'

        LE_loc = leading_edge_calc(smooth(cell_data(1,:)),xdata,0.8,0);
        IC = double(x<=LE_loc);
end


%boundary conditions
BC_x_0 = @(t) 1;
BC_x_1 = @(t) 0;

%load sparse matrices for later computation
%convection matrices
[A_pos,~,~,A_neg,~,~,D_nonlin,D_nonlin_bd,D_lin,D_lin_bd] = ...
    aMatrixConstruction(xn,dx);



tic

%run simulation
[model] = chemotaxis_computation(q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A_pos,...
    A_neg,D_nonlin,D_lin,D_nonlin_bd,D_lin_bd,x,t,xdata,tdata);
    

toc 

%plot simulation against data when t is a multiple of toCompare
toCompare = 6;
cell_data2 = cell_data(1:3*toCompare:end,1:5:end);
xdata2 = xdata(1:5:end);
tdata2 = tdata(1:3*toCompare:end);
cell_data_std = cell_data_std(1:5:end,1:3*toCompare:end);

[Xdata,Tdata] = meshgrid(xdata,tdata);
[Xint,Tint] = meshgrid(xdata2,tdata2);


udata = interp2(Xdata,Tdata,model,Xint,Tint);


figure
hold on

colors = distinguishable_colors(length(tdata2));


for i = 1:length(tdata2)
    plot(xdata2,udata(i,:),'color',colors(i,:))
    plot(xdata2,cell_data2(i,:),'.','color',[colors(i,:)],...
        'markersize',15)
    fill([xdata2 fliplr(xdata2)],[(cell_data2(i,:) +...
        cell_data_std(:,i)') fliplr(cell_data2(i,:) -...
        cell_data_std(:,i)')],colors(i,:),'Facealpha',0.1,'edgecolor','none')
end

xlabel('x')
ylabel('u(t,x)')
title(['Model ' fitpred ' to data, density ' num2str(well)])

% exportfig(gcf,[filename filename2save '.eps'],'color','rgb','fontsize',1.5)
% saveas(gcf,[filename filename2save '.fig'])

% close(gcf)

if plot_res == 1
   res_plot(udata,cell_data2,[filename filename2save]);
end

