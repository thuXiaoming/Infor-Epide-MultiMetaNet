% function [S1,I1,A1,SS,II,AA]=multiplex_meta_mc(beta1,th1,th2)
% tic
clc;clear all;close all;
N=10;T=150;time=50;
lambda1_u=0.001;ksi=0.5;mu1=0.1;
lambda2=0.05;mu2=0.15;
%lambda2=0;mu2=0;
kappa=0.5;alpha=0.3;lambda1_a=kappa*lambda1_u;
US=zeros(T+1,1);AS=zeros(T+1,1);UE=zeros(T+1,1);
AE=zeros(T+1,1);AI=zeros(T+1,1);AR=zeros(T+1,1);
n_pop=500;phi=0.5*ones(N,1);
adj=ones(N,N);
% for i=1:N*n_pop
%     adj(i,i)=0;
% end
US1=zeros(T+1,1);AS1=zeros(T+1,1);UE1=zeros(T+1,1);
AE1=zeros(T+1,1);AI1=zeros(T+1,1);AR1=zeros(T+1,1);
for tt=1:time
    tt
    A0=0.1;
    US=zeros(T+1,1);AS=zeros(T+1,1);UE=zeros(T+1,1);
    AE=zeros(T+1,1);AI=zeros(T+1,1);AR=zeros(T+1,1);
    US(1)=1-A0;AI(1)=A0;
    S_state=ones(N,n_pop);E_state=zeros(N,n_pop);
    I_state=zeros(N,n_pop);R_state=zeros(N,n_pop);
    U_state=ones(N,n_pop);A_state=zeros(N,n_pop);
    US_state=zeros(N,n_pop);AS_state=zeros(N,n_pop);
    UE_state=zeros(N,n_pop);AE_state=zeros(N,n_pop);
    AI_state=zeros(N,n_pop);AR_state=zeros(N,n_pop);
    z_vir=zeros(N,A0*n_pop);
    for i=1:N
        z_vir(i,:)=randperm(n_pop,floor(A0*n_pop)); 
    end
    for i=1:N
        for j=1:A0*n_pop
            U_state(i,z_vir(i,j))=0;A_state(i,z_vir(i,j))=1;
            S_state(i,z_vir(i,j))=0;I_state(i,z_vir(i,j))=1;
        end
    end

    for t=1:T
        %%%Information Propagation in Virtual Networks
        flag_vir_U=zeros(N,n_pop);flag_vir_A=zeros(N,n_pop);
        for i=1:N
           for j=1:n_pop
              if U_state(i,j)==1
                  flag_vir_U(i,j)=1;
                  for z=1:N
                      u1=randperm(n_pop,1);
                      %u1
                      if adj(i,z)==1&&A_state(z,u1)==1
                          if rand(1)<(lambda2*U_state(i,j))
                              flag_vir_U(i,j)=0;flag_vir_A(i,j)=1;
                              break;
                          end
                      end
                  end
%                   for z=1:N*n_pop
%                       u1=floor((z-1)/n_pop)+1;u2=z-(u1-1)*n_pop;
%                       if adj(j+(i-1)*n_pop,z)==1&&A_state(u1,u2)==1
%                           if rand(1)<(lambda2*U_state(i,j))
%                             flag_vir_U(i,j)=0;flag_vir_A(i,j)=1;
%                             break;
%                           end
%                       end
%                   end
              else
                  flag_vir_A(i,j)=1;
                  if rand(1)<mu2
                     flag_vir_U(i,j)=1;flag_vir_A(i,j)=0;
                  end
              end
           end
        end
        for i=1:N
            for j=1:n_pop
                if flag_vir_U(i,j)==1
                    U_state(i,j)=1;A_state(i,j)=0;
                else
                    U_state(i,j)=0;A_state(i,j)=1;
                end
            end
        end
        
        %%%Epidemic Spreading in Physical Networks
        flag_phy_S=zeros(N,n_pop);flag_phy_E=zeros(N,n_pop);
        flag_phy_I=zeros(N,n_pop);flag_phy_R=zeros(N,n_pop);
        loc_phy=ones(N,n_pop);

        for i=1:N
            loc_phy(i,:)=i;
        end
        %%Movement Process
        for i=1:N
           for j=1:n_pop
              if rand(1)<phi(i,1)
                  b=sample_del(N,i);
                  loc_phy(i,j)=b;
              end
           end
        end
        %%Reaction Process
        for i=1:N
            for j=1:n_pop
               if A_state(i,j)==1
                   lambda1=lambda1_a;
               else
                   lambda1=lambda1_u;
               end
               if S_state(i,j)==1
                   flag_phy_S(i,j)=1;
                   for z=1:N*n_pop
                       z1=floor((z-1)/n_pop)+1;z2=z-(z1-1)*n_pop;
                       if (loc_phy(i,j)==loc_phy(z1,z2))&&((E_state(z1,z2)==1)||(I_state(z1,z2)==1))
                          if rand(1)<(lambda1*E_state(z1,z2)*S_state(i,j)+alpha*lambda1*I_state(z1,z2)*S_state(i,j))
                             flag_phy_S(i,j)=0;flag_phy_E(i,j)=1;
                             flag_phy_I(i,j)=0;flag_phy_R(i,j)=0;
                             break;
                          end
                       end
                   end
               elseif E_state(i,j)==1
                   flag_phy_E(i,j)=1;
                   if rand(1)<ksi
                       flag_phy_S(i,j)=0;flag_phy_E(i,j)=0;
                       flag_phy_I(i,j)=1;flag_phy_R(i,j)=0;
                   end
               elseif I_state(i,j)==1
                   flag_phy_I(i,j)=1;
                   if rand(1)<mu1
                       flag_phy_S(i,j)=0;flag_phy_E(i,j)=0;
                       flag_phy_I(i,j)=0;flag_phy_R(i,j)=1;
                   end
               end
            end
        end
        
        for i=1:N
            for j=1:n_pop
                if flag_phy_S(i,j)==1
                    S_state(i,j)=1;E_state(i,j)=0;I_state(i,j)=0;R_state(i,j)=0;
                elseif flag_phy_E(i,j)==1
                    S_state(i,j)=0;E_state(i,j)=1;I_state(i,j)=0;R_state(i,j)=0;
                elseif flag_phy_I(i,j)==1
                    S_state(i,j)=0;E_state(i,j)=0;I_state(i,j)=1;R_state(i,j)=0;
                elseif flag_phy_R(i,j)==1
                    S_state(i,j)=0;E_state(i,j)=0;I_state(i,j)=0;R_state(i,j)=1;
                end
            end
        end
        %%将UI和UR状态转变为AI和AR
        for i=1:N
            for j=1:n_pop
                if U_state(i,j)==1&&(I_state(i,j)==1||R_state(i,j)==1)
                    U_state(i,j)=0;A_state(i,j)=1;
                end
            end
        end
        %%综合US-AS-UE-AE-AI-AR的状态
        for i=1:N
            for j=1:n_pop
                if U_state(i,j)==1&&S_state(i,j)==1
                    US_state(i,j)=1;AS_state(i,j)=0;UE_state(i,j)=0;
                    AE_state(i,j)=0;AI_state(i,j)=0;AR_state(i,j)=0;
                elseif A_state(i,j)==1&&S_state(i,j)==1
                    US_state(i,j)=0;AS_state(i,j)=1;UE_state(i,j)=0;
                    AE_state(i,j)=0;AI_state(i,j)=0;AR_state(i,j)=0;
                elseif U_state(i,j)==1&&E_state(i,j)==1
                    US_state(i,j)=0;AS_state(i,j)=0;UE_state(i,j)=1;
                    AE_state(i,j)=0;AI_state(i,j)=0;AR_state(i,j)=0;
                elseif A_state(i,j)==1&&E_state(i,j)==1
                    US_state(i,j)=0;AS_state(i,j)=0;UE_state(i,j)=0;
                    AE_state(i,j)=1;AI_state(i,j)=0;AR_state(i,j)=0;
                elseif A_state(i,j)==1&&I_state(i,j)==1
                    US_state(i,j)=0;AS_state(i,j)=0;UE_state(i,j)=0;
                    AE_state(i,j)=0;AI_state(i,j)=1;AR_state(i,j)=0;
                elseif A_state(i,j)==1&&R_state(i,j)==1
                    US_state(i,j)=0;AS_state(i,j)=0;UE_state(i,j)=0;
                    AE_state(i,j)=0;AI_state(i,j)=0;AR_state(i,j)=1;
                end
            end
        end
        US(t+1)=sum(US_state(:))/(N*n_pop);AS(t+1)=sum(AS_state(:))/(N*n_pop);
        UE(t+1)=sum(UE_state(:))/(N*n_pop);AE(t+1)=sum(AE_state(:))/(N*n_pop);
        AI(t+1)=sum(AI_state(:))/(N*n_pop);AR(t+1)=sum(AR_state(:))/(N*n_pop);
    end
    US1=US1+US;AS1=AS1+AS;UE1=UE1+UE;AE1=AE1+AE;AI1=AI1+AI;AR1=AR1+AR;
end
US1=US1/time;AS1=AS1/time;UE1=UE1/time;AE1=AE1/time;AI1=AI1/time;AR1=AR1/time;
ttt=0:T;
%plot(tt,p_US_ave,'-o',tt,p_AS_ave,'->',tt,p_UE_ave,'-s',tt,p_AE_ave,'-d',tt,p_AI_ave,'-*',tt,p_AR_ave,'-+')
plot(ttt,US1,'o','Linewidth',1);hold on;plot(ttt,AS1,'>','Linewidth',1);hold on;
plot(ttt,UE1,'s','Linewidth',1);hold on;plot(ttt,AE1,'d','Linewidth',1);hold on;
plot(ttt,AI1,'*','Linewidth',1);hold on;plot(ttt,AR1,'^','Linewidth',1);hold on;
        
        
%     p_US=zeros(N,T+1);p_AS=zeros(N,T+1);p_UE=zeros(N,T+1);
%     p_AE=zeros(N,T+1);p_AI=zeros(N,T+1);p_AR=zeros(N,T+1);
%     init_n_pop=500;n_pop=zeros(N,T+1);n_pop(:,1)=init_n_pop;
%     
%     I0=0.1;S_state=cell(N,T+1);E_state=cell(N,T+1);
%     I_state=cell(N,T+1);R_state=cell(N,T+1);
%     S_state{:,1}=ones(init_n_pop,1); E_state{:,1}=zeros(init_n_pop,1);
%     I_state{:,1}=zeros(init_n_pop,1);R_state{:,1}=zeros(init_n_pop,1);
%   
%     z=randperm(init_n_pop,floor(I0*init_n_pop));
%     for i=1:N
%         for j=1:floor(I0*init_n_pop)
%             S_state{i,1}(z(j))=0;E_state{i,1}(z(j))=0;
%             I_state{i,1}(z(j))=1;R_state{i,1}(z(j))=0;
%         end
%     end        
        
%         %%%Information Propagation in Virtual Networks
%         flag_vir_U=zeros(N,1);flag_vir_A=zeros(N,1);
%         for i=1:N
%             if U_state(i)==1
%                 flag_vir_U(i)=1;
%                 for j=1:N
%                     if Adj(i,j)==1&&A_state(j)==1
%                         if rand(1)<(lambda2*State(i))
%                             flag_vir_A(i)=1;flag_vir_U(i)=0;
%                             break;
%                         end
%                     end
%                 end
%             else
%                 if rand(1)<mu2
%                     flag_vir_A(i)=0;flag_vir_U(i)=1;
%                 end
%             end
%         end
%         for i=1:N
%             if flag_vir_U(i)==1
%                 U_state(i)=1;A_state(i)=0;
%             else
%                 U_state(i)=0;A_state(i)=1;
%             end
%         end
%         U(t+1)=sum(U_state)/N;A(t+1)=sum(A_state)/N;
            

