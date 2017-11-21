%%%% chemotaxis_simulate written 11-7-17 by JTN 

% define x, t grids
xn = 101;
dt = 1e-4;

x = linspace(0,1,xn);
t = 0:dt:.5;

dx = x(2) - x(1);
tn = length(t);


%params

D = .05;
Dv = 0.1;
chi = 0.05;
alpha = 0;
beta = 2;

q = [D,chi,alpha,beta,Dv];

%interior, boundary points
[x_int,xbd_0,xbd_1] = int_bd_def(xn);

ICu = double(x<=.5);
ICv = 1*ones(size(x));%+1/10*exp(-10*(x).^2);


IC = [ICu(:) ; ICv(:)];

[A_pos,~,~,A_neg,~,~,D_nonlin,D_nonlin_bd,D_lin,...
    D_lin_bd] = aMatrixConstruction(xn,q,dx);

tic

    model = chemotaxis_computation(q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A_pos,...
        A_neg,D_nonlin,D_lin,D_nonlin_bd,D_lin_bd);

toc

u = model(1:xn,:);
v = model(xn+1:end,:);

count = 1;
for i = 1:floor(tn/50):tn
    hold off
    plot(x,u(:,i))
    hold on
    plot(x,v(:,i))
    axis([0 1 0 1.2])
    title(['t = ' num2str(t(i))])
    xlabel('x')
%     ylabel('u')
    legend('u','v')
    pause(.5)
%     exportfig(gcf,['chem_fig_' num2str(count) '.eps'],'color','rgb','fontsize',1.5)
%     count = count + 1;

end