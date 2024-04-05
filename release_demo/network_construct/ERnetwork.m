function [m,Adj,X,Y]=erdosRenyi(nv,p)
center=[0,0];
theta=linspace(0,2*pi,nv);%devide the circle into 'nv' parts to put the vertices onto it! using polar coordinate!
rho=ones(1,nv);%the more nv, the larger the circle radius!
[X,Y] = pol2cart(theta',rho');%transfer polar cooedinate to decart coordinate!
X=X+center(1);
Y=Y+center(2);%Registering the TWO new vectors coordinates by X,Y!
P=nchoosek(1:nv,2);%choose two elements from 1 to nv, registered as a matrix of nv*2;
C=(rand(length(P),1)<=p);% Randomly generate a P_equal-length column vector and reserving the necessary links should added
                        % which are determinied by parability 'p', then registered it by 'C'.
Adj=sparse(P(:,1),P(:,2),C);% Creat the Adjacent Materix by a sparse matrix type.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the Adj to a Symmetry Square Matrix! 2005-04-19
%Symmetry a non-square Matrix.
m=full(Adj);
m(nv,nv)=[0];
for i=1:nv
   for j=1:nv
       if (i>=j)
          m(j,i)=m(i,j);
          m(i,j)=m(i,j);
       else
       end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting section (comment to avoid plotting)
f=figure;
set(f,'Color','w')
gplot(Adj,[X Y],'k-')
hold on;
h=scatter(X,Y,5,'y','filled');
set(gca,'Visible','Off');
set(h,'LineWidth',2,'MarkerEdgeColor','k');
axis square
%end of plotting section