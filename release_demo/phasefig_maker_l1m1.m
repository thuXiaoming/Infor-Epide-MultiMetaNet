%%%%本程序实现的是lambda1_mu1的相图绘制，包含阈值曲线
clc;clear all;close all;

load('D:\Code_multiplex_metapop\lambda_c.mat')
phi_0=0.5;i=0;j=0;N1=50;N=10;
p_US_final=zeros(N1,N1);p_AS_final=zeros(N1,N1);p_UE_final=zeros(N1,N1);
p_AE_final=zeros(N1,N1);p_AI_final=zeros(N1,N1);p_AR_final=zeros(N1,N1);
X1=0.0001:0.0002:(0.0001+0.0002*(N1-1));Y1=(0.05+0.01*(N1-1)):-0.01:0.05;
for mu1=(0.05+0.01*(N1-1)):-0.01:0.05
    i=i+1;j=0;
    for lambda1_u=0.0001:0.0002:(0.0001+0.0002*(N1-1))
        j=j+1;
        [p_US_final(i,j),p_AS_final(i,j),p_UE_final(i,j),p_AE_final(i,j),p_AI_final(i,j),p_AR_final(i,j)]=emu_para(lambda1_u,mu1,phi_0);
    end
end
xmax=max(X1);xmin=min(X1);
ymax=max(Y1);ymin=min(Y1);
[X,Y]=meshgrid(linspace(xmin,xmax,N1),linspace(ymin,ymax,N1));
Z1 = griddata(X1,Y1,p_AR_final,X,Y,'v4');
% Z2 = griddata(X1,Y1,p_US_final,X,Y,'v4');
% Z3 = griddata(X1,Y1,p_AS_final,X,Y,'v4');
% Z4 = griddata(X1,Y1,p_UE_final,X,Y,'v4');
% Z5 = griddata(X1,Y1,p_AE_final,X,Y,'v4');
% Z6 = griddata(X1,Y1,p_AI_final,X,Y,'v4');

%%% 等高线法
figure('NumberTitle','off','Name','等高线法_AR','Color','w');
contourf(X,Y,Z1,N1, 'LineColor','none');
colormap('jet');
colorbar;
hold on;
mu_mid=0.05:0.01:(0.05+0.01*(N1-1));
plot(lambda_c,mu_mid,'w-.x','linewidth',1.5);
