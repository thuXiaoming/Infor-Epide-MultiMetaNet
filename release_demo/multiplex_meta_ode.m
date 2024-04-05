clc;clear all;close all;

N=10;   %total number of patches(cities)
T=150;   %dynamic evolution period
p_US=zeros(N,T+1);p_AS=zeros(N,T+1);p_UE=zeros(N,T+1);
p_AE=zeros(N,T+1);p_AI=zeros(N,T+1);p_AR=zeros(N,T+1);
n_pop=zeros(N,T+1); %population size of each patch
p_US(:,1)=0.9;p_AS(:,1)=0;p_UE(:,1)=0;
p_AE(:,1)=0;p_AI(:,1)=0.1;p_AR(:,1)=0;
n_pop(:,1)=500;epsilon_0=0.3;

lambda1_u=0.001;ksi=0.5;mu1=0.1;   %spread parameters of pyhsical layers 
lambda2=0.05;mu2=0.15;           %spread parameters of virtual layers
%lambda2=0;mu2=0;

kappa=0.5;alpha=0.3;
lambda1_a=kappa*lambda1_u;
phi=0.5*ones(N,1);  %mobility rate of each patch
Phi_U=zeros(N,T);Phi_A=zeros(N,T);n_jtoi=zeros(N,N,T);
q_U=zeros(N,T);q_A=zeros(N,T);g=zeros(N,T);
gamma=ones(N,1);
epsilon=zeros(N,N);
R=ones(N,N)/(N-1);w_v=ones(N,N);
for i=1:N
    R(i,i)=0;
    %w_v(i,i)=0;
end

for t=1:T
    t
    for i=1:N
        pai_U=1;pai_A=1;pai_g=1;
        for j=1:N
            %n_jtoi(j,i,t)=(KronDel(i,j)*(1-phi(j,1))*n_pop(j,t)+phi(j,1)*R(j,i)*n_pop(j,t))/n_pop(j,t);
            n_jtoi(j,i,t)=KronDel(i,j)*(1-phi(j,1))*n_pop(j,t)+phi(j,1)*R(j,i)*n_pop(j,t);
            pai_U=pai_U*((1-lambda1_u*(p_UE(j,t)+p_AE(j,t))-alpha*lambda1_u*p_AI(j,t))^n_jtoi(j,i,t));
            pai_A=pai_A*((1-lambda1_a*(p_UE(j,t)+p_AE(j,t))-alpha*lambda1_a*p_AI(j,t))^n_jtoi(j,i,t));
            pai_g=pai_g*(1-lambda2*w_v(i,j)*(p_AS(j,t)+p_AE(j,t)+p_AI(j,t)+p_AR(j,t)));
        end
        Phi_U(i,t)=1-pai_U;Phi_A(i,t)=1-pai_A;g(i,t)=1-pai_g;
    end
    for i=1:N
        B_mid_U=0;B_mid_A=0;
        for j=1:N
            B_mid_U=B_mid_U+R(i,j)*Phi_U(j,t);
            B_mid_A=B_mid_A+R(i,j)*Phi_A(j,t);
        end
        q_U(i,t)=(1-phi(i,1))*Phi_U(i,t)+phi(i,1)*B_mid_U;
        q_A(i,t)=(1-phi(i,1))*Phi_A(i,t)+phi(i,1)*B_mid_A;
    end
    %%dynamics of population size in each patch
    for i=1:N
        C_mid1=0;C_mid2=0;
        for j=1:N
            if (p_AI(i,t)^gamma(i,1)-p_AI(j,t)^gamma(j,1))>exp(-10)
                epsilon(i,j)=epsilon_0;
            else
                epsilon(i,j)=0;
            end
            C_mid1=C_mid1+R(i,j)*epsilon(i,j);C_mid2=C_mid2+R(j,i)*epsilon(j,i)*phi(j,1)*n_pop(j,t);
        end
        n_pop(i,t+1)=n_pop(i,t)-phi(i,1)*n_pop(i,t)*C_mid1+C_mid2;
    end
    
    %%Spreading Process Evolution Iteration
    for i=1:N
        p_US(i,t+1)=p_US(i,t)*(1-g(i,t))*(1-q_U(i,t))+p_AS(i,t)*mu2*(1-q_U(i,t));
        p_AS(i,t+1)=p_US(i,t)*g(i,t)*(1-q_A(i,t))+p_AS(i,t)*(1-mu2)*(1-q_A(i,t));
        p_UE(i,t+1)=p_US(i,t)*(1-g(i,t))*q_U(i,t)+p_AS(i,t)*mu2*q_U(i,t)+p_UE(i,t)*(1-g(i,t))*(1-ksi)+p_AE(i,t)*mu2*(1-ksi);
        p_AE(i,t+1)=p_US(i,t)*g(i,t)*q_A(i,t)+p_AS(i,t)*(1-mu2)*q_A(i,t)+p_UE(i,t)*g(i,t)*(1-ksi)+p_AE(i,t)*(1-mu2)*(1-ksi);
        p_AI(i,t+1)=(p_UE(i,t)+p_AE(i,t))*ksi+p_AI(i,t)*(1-mu1);
        p_AR(i,t+1)=p_AR(i,t)+p_AI(i,t)*mu1;
    end
    
end
AA=p_US+p_AS+p_UE+p_AE+p_AI+p_AR; 
p_US_ave=sum(p_US,1)/N;p_AS_ave=sum(p_AS,1)/N;p_UE_ave=sum(p_UE,1)/N;
p_AE_ave=sum(p_AE,1)/N;p_AI_ave=sum(p_AI,1)/N;p_AR_ave=sum(p_AR,1)/N;
tt=0:T;
%plot(tt,p_US_ave,'-o',tt,p_AS_ave,'->',tt,p_UE_ave,'-s',tt,p_AE_ave,'-d',tt,p_AI_ave,'-*',tt,p_AR_ave,'-+')
plot(tt,p_US_ave,'-','Linewidth',1);hold on;plot(tt,p_AS_ave,'-','Linewidth',1);hold on;
plot(tt,p_UE_ave,'-','Linewidth',1);hold on;plot(tt,p_AE_ave,'-','Linewidth',1);hold on;
plot(tt,p_AI_ave,'-','Linewidth',1);hold on;plot(tt,p_AR_ave,'-','Linewidth',1);hold on;
