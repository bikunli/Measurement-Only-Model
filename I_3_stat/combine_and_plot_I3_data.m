clear 
close all

MutualInfo3_v_L32 = [];
load('compute_I3_StartingPure_v1p0_L32.mat');
MutualInfo3_v_L32 = [MutualInfo3_v_L32;MutualInfo3_mat];
load('compute_I3_StartingPure_v1p1_L32.mat');
MutualInfo3_v_L32 = [MutualInfo3_v_L32;MutualInfo3_mat];
MutualInfo3_v_L32 = mean(MutualInfo3_v_L32,1);

qx_v_L32 = qx_v;

%---------
MutualInfo3_v_L64 = [];
load('compute_I3_StartingPure_v1p0_L64.mat');
MutualInfo3_v_L64 = [MutualInfo3_v_L64;MutualInfo3_mat];
load('compute_I3_StartingPure_v1p1_L64.mat');
MutualInfo3_v_L64 = [MutualInfo3_v_L64;MutualInfo3_mat];
load('compute_I3_StartingPure_v1p2_L64.mat');
MutualInfo3_v_L64 = [MutualInfo3_v_L64;MutualInfo3_mat];
MutualInfo3_v_L64 = mean(MutualInfo3_v_L64,1);

qx_v_L64 = qx_v;

%---------

MutualInfo3_v_L128 = [] ;
load('compute_I3_StartingPure_v1p0_L128.mat');
MutualInfo3_v_L128 = [MutualInfo3_v_L128;MutualInfo3_mat];
load('compute_I3_StartingPure_v1p1_L128.mat');
MutualInfo3_v_L128 = [MutualInfo3_v_L128;MutualInfo3_mat];
for i = 0:4
    load(['compute_I3_StartingPure_v2p',num2str(i),'_L128.mat']);
    MutualInfo3_v_L128 = [MutualInfo3_v_L128;MutualInfo3_mat];
end

MutualInfo3_v_L128 = mean(MutualInfo3_v_L128,1);

qx_v_L128 = qx_v;

%---------

MutualInfo3_v_L256 = [] ;
for i = 0:60
    load(['compute_I3_StartingPure_v1p',num2str(i),'_L256.mat']);
    MutualInfo3_v_L256 = [MutualInfo3_v_L256;MutualInfo3_mat];
end

MutualInfo3_v_L256 = mean(MutualInfo3_v_L256,1);

qx_v_L256 = qx_v;

%%%========================================================
DC = DefaultColor();

subplot(211)
plot(qx_v_L32, MutualInfo3_v_L32,'-*','MarkerSize',4,'Color',0.65*[1,1,1],'DisplayName','$L = 32$'); hold on
plot(qx_v_L64, MutualInfo3_v_L64,'-d','MarkerSize',4,'Color',DC(1,:),'DisplayName','$L = 64$'); hold on
plot(qx_v_L128, MutualInfo3_v_L128,'-v','MarkerSize',4,'Color',DC(2,:),'DisplayName','$L = 128$'); hold on
plot(qx_v_L256, MutualInfo3_v_L256,'-o','MarkerSize',4,'Color',DC(5,:),'DisplayName','$L = 256$'); hold on

legend('interpreter','latex','FontSize',12,'Location','SouthWest');
grid on
xlabel('$q_X$','Interpreter','latex','FontSize',16);
ylabel('$\overline{\mathcal{I}_3}$','Interpreter','latex','FontSize',16);
xlim([0.18,0.35]);
ylim([-4,2.25]);



subplot(212)
nu = 1.1;
plot((qx_v_L32 - 0.274)*32^(1/nu), MutualInfo3_v_L32,'-*','MarkerSize',4,'Color',0.65*[1,1,1],'DisplayName','$L = 32$'); hold on
plot((qx_v_L64 - 0.274)*64^(1/nu), MutualInfo3_v_L64,'-d','MarkerSize',4,'Color',DC(1,:),'DisplayName','$L = 64$'); hold on
plot((qx_v_L128 - 0.274)*128^(1/nu), MutualInfo3_v_L128,'-v','MarkerSize',4,'Color',DC(2,:),'DisplayName','$L = 128$'); hold on
plot((qx_v_L256 - 0.274)*256^(1/nu), MutualInfo3_v_L256,'-o','MarkerSize',4,'Color',DC(5,:),'DisplayName','$L = 256$'); hold on

legend('interpreter','latex','FontSize',12,'Location','SouthWest');
grid on
xlabel('$(q_X - q_{X,c})L^{1/\nu}$','Interpreter','latex','FontSize',16);
ylabel('$\overline{\mathcal{I}_3}$','Interpreter','latex','FontSize',16);
xlim([-3,2]);
ylim([-2,2]);


function [dc] = DefaultColor()

dc = [         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];


end