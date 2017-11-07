%%%% chemotaxis_simulate written 11-7-17 by JTN 

% define x, t grids
xn = 201;
dt = 1e-3;

x = linspace(0,1,xn);
t = 0:dt:10;

dx = x(2) - x(1);
tn = length(t);


%params

D = 0;
chi = 0.05;
n = 0;
alpha = 1;

q = [D,chi,n,alpha];

%interior, boundary points
[x_int,xbd_0,xbd_1] = int_bd_def(xn);

ICu = ones(xn,1);
ICv = 1+1/10*exp(-10*(x-1).^2);

IC = [ICu(:) ; ICv(:)];

[A_pos,~,~,A_neg,~,~,D_nonlin,D_nonlin_bd,D_lin,...
    D_lin_bd] = aMatrixConstruction(xn,D,n,dx);

model = chemotaxis_computation(q,dx,xn,x_int,xbd_0,xbd_1,dt,tn,IC,A_pos,...
    A_neg,D_nonlin,D_lin,D_nonlin_bd,D_lin_bd);

u = model(1:xn,:);
v = model(xn+1:end,:);

for i = 1:10
    hold off
    plot(u(:,i))
    hold on
    plot(v(:,i))
    axis([0 201 .8 1.2])
    pause(.125)
end