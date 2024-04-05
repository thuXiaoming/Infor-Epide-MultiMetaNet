clc;
clear all;

m0=input('请输入未增长前的网络节点个数m0：');
m=input('请输入每次引入新节点时新生成的边数m：');
N=input('请输入增长后的网络节点总数N：');
disp('初始网络时m0个节点的连接情况：1表示都是孤立点；2表示构成完全图；3表示随机连接一些边');
se=input('请选择初始网络情况1，2或3:');
if m>m0
    disp('输入参数m不合法');return;
end
x=100*rand(1,m0);y=100*rand(1,m0);      %构造初始用于画图的m0个节点坐标
if se==1
    A=zeros(m0);
elseif se==2
    A=ones(m0);A(1:m0+1:m0^2)=0;        %对角线元素置0
else
    A=zeros(m0);B=rand(m0);B=tril(B);   %截取下三角元素
    A(B<=0.4 & B>0)=1;            %按照0.1概率进行连边
    A=A+A';         %构造完整的邻接矩阵
end
for k=m0+1:N
    x(k)=100*rand;y(k)=100*rand;        %生成用于当前节点画图的坐标
    p=(sum(A)+1)/sum(sum(A)+1);         %计算所有节点的连接概率
    pp=cumsum(p);           %求累计分布
    A(k,k)=0;               %加入新的连边之前，邻接矩阵扩充维数
    ind=[];                 %新节点所连节点的初始集合
    while length(ind)<m
        jj=find(pp>rand);jj=jj(1);  %用赌轮法选择连边节点的编号
        ind=union(ind,jj);          %使用union保证选择的节点不重复
    end
    A(k,ind)=1;A(ind,k)=1;          %构造加边以后新的邻接矩阵
end
for n=1:N
    A(n,n)=0;
end
plot(x,y,'ro','MarkerEdgeColor','g','MarkerFaceColor','r','markersize',8);
hold on,A2=tril(A);[i,j]=find(A2);      %找邻接矩阵下三角元素的非零元素
for k=1:length(i)
    plot([x(i(k)),x(j(k))],[y(i(k)),y(j(k))],'linewidth',1.2);
end
%mydegree(A);

    

    