function verify(appcase,para)

figure

switch appcase
    case  'fp64'
        h.X = load('./data/admm/admm_xk_fp64.mat')';
        h.X = h.X.x64;
        h.Z = load('./data/admm/admm_zk_fp64.mat')';
        h.Z = h.Z.z64;
        h.U = load('./data/admm/admm_uk_fp64.mat')';
        h.U = h.U.u64;
        xadmm_optval = load('./data/admm/admm_fk_fp64.mat')';
        xadmm_optval = xadmm_optval.y64;
        % ground truth for fp64
        ground_x = load('./data/admm/admm_xk_ref_fp64.mat')';
        ground_x = ground_x.x64_ref;
        ground_z = load('./data/admm/admm_zk_ref_fp64.mat')';
        ground_z = ground_z.z64_ref;
        ground_u = load('./data/admm/admm_uk_ref_fp64.mat')';
        ground_u = ground_u.u64_ref;
        admm_optval = load('./data/admm/admm_fk_ref_fp64.mat')';
        admm_optval = admm_optval.y64_ref;
    case  'fp32'
        h.X = load('./data/admm/admm_xk_fp32.mat')';
        h.X = h.X.x32;
        h.Z = load('./data/admm/admm_zk_fp32.mat')';
        h.Z = h.Z.z32;
        h.U = load('./data/admm/admm_uk_fp32.mat')';
        h.U = h.U.u32;
        xadmm_optval = load('./data/admm/admm_fk_fp32.mat')';
        xadmm_optval = xadmm_optval.y32;
        % ground truth for fp32
        ground_x = load('./data/admm/admm_xk_ref_fp32.mat')';
        ground_x = ground_x.x32_ref;
        ground_z = load('./data/admm/admm_zk_ref_fp32.mat')';
        ground_z = ground_z.z32_ref;
        ground_u = load('./data/admm/admm_uk_ref_fp32.mat')';
        ground_u = ground_u.u32_ref;
        admm_optval = load('./data/admm/admm_fk_ref_fp32.mat')';
        admm_optval = admm_optval.y32_ref;
    
    case  'fx16'
        h.X = load('./data/admm/admm_xk_fx16.mat')';
        h.X = h.X.x16;
        h.Z = load('./data/admm/admm_zk_fx16.mat')';
        h.Z = h.Z.z16;
        h.U = load('./data/admm/admm_uk_fx16.mat')';
        h.U = h.U.u16;
        xadmm_optval = load('./data/admm/admm_fk_fx16.mat')';
        xadmm_optval = xadmm_optval.y16;
        % ground truth for fx16
        ground_x = load('./data/admm/admm_xk_ref_fx16.mat')';
        ground_x = ground_x.x16_ref;
        ground_z = load('./data/admm/admm_zk_ref_fx16.mat')';
        ground_z = ground_z.z16_ref;
        ground_u = load('./data/admm/admm_uk_ref_fx16.mat')';
        ground_u = ground_u.u16_ref;
        admm_optval = load('./data/admm/admm_fk_ref_fx16.mat')';
        admm_optval = admm_optval.y16_ref;
    case  'fx12'
        h.X = load('./data/admm/admm_xk_fx12.mat')';
        h.X = h.X.x12;
        h.Z = load('./data/admm/admm_zk_fx12.mat')';
        h.Z = h.Z.z12;
        h.U = load('./data/admm/admm_uk_fx12.mat')';
        h.U = h.U.u12;
        xadmm_optval = load('./data/admm/admm_fk_fx12.mat')';
        xadmm_optval = xadmm_optval.y12;
        % ground truth for fx12
        ground_x = load('./data/admm/admm_xk_ref_fx12.mat')';
        ground_x = ground_x.x12_ref;
        ground_z = load('./data/admm/admm_zk_ref_fx12.mat')';
        ground_z = ground_z.z12_ref;
        ground_u = load('./data/admm/admm_uk_ref_fx12.mat')';
        ground_u = ground_u.u12_ref;
        admm_optval = load('./data/admm/admm_fk_ref_fx12.mat')';
        admm_optval = admm_optval.y12_ref;
    case  'fx8'
        h.X = load('./data/admm/admm_xk_fx8.mat')';
        h.X = h.X.x8;
        h.Z = load('./data/admm/admm_zk_fx8.mat')';
        h.Z = h.Z.z8;
        h.U = load('./data/admm/admm_uk_fx8.mat')';
        h.U = h.U.u8;
        xadmm_optval = load('./data/admm/admm_fk_fx8.mat')';
        xadmm_optval = xadmm_optval.y8;
        % ground truth for fx8
        ground_x = load('./data/admm/admm_xk_ref_fx8.mat')';
        ground_x = ground_x.x8_ref;
        ground_z = load('./data/admm/admm_zk_ref_fx8.mat')';
        ground_z = ground_z.z8_ref;
        ground_u = load('./data/admm/admm_uk_ref_fx8.mat')';
        ground_u = ground_u.u8_ref;
        admm_optval = load('./data/admm/admm_fk_ref_fx8.mat')';
        admm_optval = admm_optval.y8_ref;
    case  'fx2'
        h.X = load('./data/admm/admm_xk_fx2.mat')';
        h.X = h.X.x2;
        h.Z = load('./data/admm/admm_zk_fx2.mat')';
        h.Z = h.Z.z2;
        h.U = load('./data/admm/admm_uk_fx2.mat')';
        h.U = h.U.u2;
        xadmm_optval = load('./data/admm/admm_fk_fx2.mat')';
        xadmm_optval = xadmm_optval.y2;
        % ground truth for fx2
        ground_x = load('./data/admm/admm_xk_ref_fx2.mat')';
        ground_x = ground_x.x2_ref;
        ground_z = load('./data/admm/admm_zk_ref_fx2.mat')';
        ground_z = ground_z.z2_ref;
        ground_u = load('./data/admm/admm_uk_ref_fx2.mat')';
        ground_u = ground_u.u2_ref;
        admm_optval = load('./data/admm/admm_fk_ref_fx2.mat')';
        admm_optval = admm_optval.y2_ref;
    case  'fpx'
        h.X = load('./data/admm/admm_xk_fpx.mat')';
        h.Z = load('./data/admm/admm_zk_fpx.mat')';
        h.U = load('./data/admm/admm_uk_fpx.mat')';
    case  'fxp'
        h.X = load('./data/admm/admm_xk_fxp.mat')';
        h.Z = load('./data/admm/admm_zk_fxp.mat')';
        h.U = load('./data/admm/admm_uk_fxp.mat')';
    case  'spx1'
        h.X = load('./data/admm/admm_xk_spx1.mat')';
        h.Z = load('./data/admm/admm_zk_spx1.mat')';
        h.U = load('./data/admm/admm_uk_spx1.mat')';
    case  'spx2'
        h.X = load('./data/admm/admm_xk_spx2.mat')';
        h.Z = load('./data/admm/admm_zk_spx2.mat')';
        h.U = load('./data/admm/admm_uk_spx2.mat')';
end
        
% [para] = loaddata_admm(para);
% 
% admm_optval = objective_lasso(para.A, para.yq(end,:)', para.lambda_yq_dct, ground_x(:,end), ground_z(:,end));
% 
% for i=1:size(h.X,2)
%     xadmm_optval = objective_lasso(para. A, para.yq(end,:)', para.lambda_yq_dct, h.X(:,i), h.Z(:,i));
% end

n=size(h.X,1);
n2 = min(size(h.X,2),size(ground_x,2));

ground_x = ground_x(:,1:n2);
ground_z = ground_z(:,1:n2);
h.X = h.X(:,1:n2);
h.Z = h.Z(:,1:n2);

% x0=zeros(n,1);
% z0=zeros(n,1);
% u0=zeros(n,1);

x0=h.X(:,1);
z0=h.Z(:,1);
u0=zeros(n,1);

% x_opt=ground_x(:,end);
% z_opt=ground_z(:,end);
x_opt=h.X(:,end);
z_opt=h.Z(:,end);
A=eye(n);
B=-eye(n);
c=0;

% rho=para.rho;
probLambda= para.gamma;
paramLambda= 1/para.rho;

% tolerances
%h.SolvTol1=max(h.X-ground_x);
h.SolvTol1=vecnorm(h.X-ground_x,2).^2/2;
%h.SolvTol2=max(h.Z-ground_z);
h.SolvTol2=vecnorm(h.Z-ground_z,2).^2/2;
admm_bound
end