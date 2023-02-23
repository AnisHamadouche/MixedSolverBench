% x0=sprandn(n,1,0.1);
% z0=sprandn(n,1,0.1);
% u0=zeros(n,1);
% x_opt=x_opt;
% 
% knorm=20;
% A=eye(n);
% B=-eye(n);
% c=0;
% rho=1.2;
% probLambda=1;
% paramLambda=1/rho;
% skip_r=50;
% skip_l=50;
% 
% h_ground = solv_ks_admm(x0, mA, mb, A,B,c,probLambda,knorm,x0,z0,u0,'WL',eye(n),'rho',rho,'skip_r',1,'skip_l',1);
% h = solv_ks_admm(x0, mA, mb, A,B,c,probLambda,knorm,x0,z0,u0,'WL',eye(n),'rho',rho,'skip_r',skip_r,'skip_l',skip_l);


zeta = 20*sqrt(log(2)); %probability parameter
% admm_optval = h_ground.admm_optval;
% xadmm_optval = h.admm_optval;
X_iter = h.X;
Z_iter = h.Z;
U_iter = h.U;
MAX_ITER = size(X_iter,2);
% x0 = X_iter(:,1);
% z0 = Z_iter(:,1);
% u0 = U_iter(:,1);
%SlvTol1=repmat(sqrt(2*h.SolvTol1/(beta)),1,MAX_ITER);
SlvTol1=h.SolvTol1;
%SlvTol2=repmat(uperrbound(knorm,beta,skip_r,skip_l,n),1,MAX_ITER);
SlvTol2=h.SolvTol2;
eps_g = SlvTol1; %tolerance of the solution of x subproblem
eps_h = SlvTol2; %tolerance of the solution of z subproblem

Sigma = eye(n)-(1-1/paramLambda)*eye(n);
paramL = eye(n);
eig_max = max(eig(Sigma'*Sigma));
epsilon_0 = max(max(SlvTol1),max(SlvTol2)); %max()
%epsilon_0 = max(median(SlvTol1),median(SlvTol2)); %max()
k = 1;
prob_bound = zeros(1,MAX_ITER);
prob_bound(k) = ((1-1/paramLambda)*norm(x0-x_opt)^2+norm(z0-z_opt)^2)/(2);%+dot(paramL*(U_iter(:,2)+B*(Z_iter(:,2)-Z_iter(:,1))),U_iter(:,2)-U_iter(:,1))/(2)/paramLambda;...%(2*(k+1))
    %+(sqrt(2*eps_g(k))*norm(Z_iter(:,k)-z_opt))/(k+1)...
    %+(sqrt(2*eps_h(k))*norm(X_iter(:,k)-x_opt))/(k+1)...%+(sqrt(2*eig_max*eps_h(k))*norm(X(:,k)-x_opt))/(k+1)...
    %-(zeta)*epsilon_0/sqrt(k+1);
for k=2:MAX_ITER
    prob_bound(k) = mean(eps_g)+mean(eps_h)...%(sum(eps_g(1:k))+sum(eps_h(1:k)))/(k+1)...
    +((1-1/paramLambda)*norm(x0-x_opt)^2+norm(z0-z_opt)^2)/(2*(k+1))...
    +(sqrt(2*eps_h(k))*norm(Z_iter(:,k)-z_opt))/(k+1)...
    +(sqrt(2*eig_max*eps_g(k))*norm(X_iter(:,k)-x_opt))/(k+1)...%+(sqrt(2*eps_g(k))*norm(X_iter(:,k)-x_opt))/(k+1)...
    +(zeta)*epsilon_0/sqrt(k+1)...
    -dot(paramL*(U_iter(:,k)+B*(Z_iter(:,k)-Z_iter(:,k-1))),U_iter(:,k)-U_iter(:,k-1))/(k+1)/paramLambda;
    %+vecnorm(paramL*U_iter(:,k)+B*(Z_iter(:,k)-Z_iter(:,k-1)),2)*vecnorm(U_iter(:,k)-U_iter(:,k-1),2)/(k+1)/paramLambda;
    %+abs(dot(paramL*U_iter(:,k)+B*(Z_iter(:,k)-Z_iter(:,k-1)),U_iter(:,k)-U_iter(:,k-1))/(k+1)/paramLambda);
    %+vecnorm(paramL*U_iter(:,k)+B*(Z_iter(:,k)-Z_iter(:,k-1)),2)*vecnorm(U_iter(:,k)-U_iter(:,k-1),2)/(k+1)/paramLambda;    
end

eps_g = zeros(1,MAX_ITER);
eps_h = zeros(1,MAX_ITER);
epsilon_0 = 0;
k = 1;
bound0 = zeros(1,MAX_ITER);
bound0(k) = ((1-1/paramLambda)*norm(x0-x_opt)^2+norm(z0-z_opt)^2)/(2);%+dot(paramL*(U_iter(:,2)+B*(Z_iter(:,2)-Z_iter(:,1))),U_iter(:,2)-U_iter(:,1))/(2)/paramLambda;%(2*(k+1));
for k=2:MAX_ITER
    bound0(k) = ((1-1/paramLambda)*norm(x0-x_opt)^2+norm(z0-z_opt)^2)/(2*(k+1))...%+vecnorm(paramL*U_iter(:,k)+Z_iter(:,k)-Z_iter(:,k-1),2)*vecnorm(U_iter(:,k)-U_iter(:,k-1),2)/(k+1)/paramLambda;
    -dot(paramL*U_iter(:,k)+Z_iter(:,k)-Z_iter(:,k-1),U_iter(:,k)-U_iter(:,k-1))/(k+1)/paramLambda;
end

%Calculate Probability
% CMP=[xadmm_optval(1:end)-xadmm_optval(end)]>prob_bound;
% CMP=sum(CMP(:)==1)/size(CMP,2);
% fprintf('Theoretical Probability Rate =%.2f\nExperimental Rate =%.2f\n', 1-2*exp(-zeta^2/2), 1-CMP)

%Plot
%semilogy(abs(xadmm_optval(1:end)-admm_optval(end)),'Color',[0.80,0.80,0.80],'linewidth',2,'DisplayName','$f^k - f^\star$')
%loglog(abs(xadmm_optval(1:end)-admm_optval(end)),'Color',[0.80,0.80,0.80],'linewidth',2,'DisplayName','$f^k - f^\star$')
%loglog(abs(xadmm_optval(1:end)-admm_optval(end)),'Color',[0.00,0.45,0.74],'linewidth',3,'DisplayName','$f^k - f^\star$')
xadmm_avg = [];
for k=1:1000
    xadmm_avg=[xadmm_avg, sum(xadmm_optval(1:k))/(k+1)];
end
loglog(abs(xadmm_avg(1:end)-admm_optval(end)),'Color',[0.00,0.45,0.74],'linewidth',3,'DisplayName','$f^k - f^\star$')
hold
% semilogy(bound0(1:end),'r','linewidth',0.5,'DisplayName','Error free')
% semilogy(prob_bound(1:end),':', 'linewidth',3,'Color',[0 0 0],'DisplayName','Theorem 3 ($\gamma = 2\sqrt{log(2)}$)')
loglog(bound0(1:end),'r','linewidth',0.5,'DisplayName','Error free')
loglog(prob_bound(1:end),'-', 'linewidth',3,'Color',[0 0 0],'DisplayName','Theorem 3 ($\gamma = 20\sqrt{log(2)}$)')
ylabel('$f^k - f^\star$','interpreter','latex');
xlabel('$Iterations, k$','interpreter','latex');
xlim([1 k])
hl = legend('show');
set(hl, 'Interpreter','latex')
grid on

mean_end = mean(xadmm_optval(end)-admm_optval(end));
subopt_pred = abs(double((prob_bound(end)-mean_end)))
subopt_pred0 = abs(double((bound0(end)-mean_end)))

%Plot normal
% plot(xadmm_optval(2:end)-admm_optval(end),'b','linewidth',2,'DisplayName','$f^k - f^\star$')
% hold
% plot(bound0(2:end),'g','linewidth',0.5,'DisplayName','Error free')
% plot(prob_bound(2:end),':', 'linewidth',3,'Color',[0 0 0],'DisplayName','Theorem 3')
% ylabel('$f^k - f^\star$','interpreter','latex');
% xlabel('$Iterations, k$','interpreter','latex');
% xlim([1 k])
% hl = legend('show');
% set(hl, 'Interpreter','latex')
% grid on