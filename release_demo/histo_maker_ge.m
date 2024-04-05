%%%%本程序用来实现每个patch的异构的phi vs. gamma
clc;clear all;close all;
N=10;
p_AR_final=zeros(5,N);p_AR_ave_final=zeros(5,1);

for mode=1:5
    [p_AR_final(mode,:),p_AR_ave_final(mode,1)]=histo_ge(mode);
end

pp=[p_AR_ave_final,p_AR_final];
figureHandle = figure;
b = bar3(pp,0.75);
% 赋色
for k = 1:length(b)
   zdata = b(k).ZData;
   b(k).CData = zdata;
   b(k).FaceColor = 'interp';
end
set(gca,'Xticklabel',{'Ave','1','2','3','4','5','6','7','8','9','10'})
xlabel('patch number')
set(gca,'Yticklabel',{'1','2','3','4','5'})
ylabel('mode')
zlabel('recovered fraction')

colormap('jet')