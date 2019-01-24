close all
clear all
%close all
s1='bin/coloring_';
s3='.txt';
for jj=0:1
s2=num2str(jj);
s=strcat(s1,s2,s3)
A=load(s);   
n=length(A(:,1));
radius=1/(2*n);
local_color{1}=[0,0,0];
local_color{2}=[0,1,0];
local_color{3}=[0,0,1];
local_color{4}=[1,1,0.25];
local_color{5}=[1.,0.0,0.75];
 figure
hold on

for ii=1:n
center=[A(ii,1),A(ii,2)];
% if(jj==0)
txt = num2str(A(ii,4));
text(center(1),center(2),txt)
        if(A(ii,3)==0)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','y','LineWidth',2);
        elseif(A(ii,3)==1)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','g','LineWidth',2);
        elseif(A(ii,3)==2)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','b','LineWidth',2);
        elseif(A(ii,3)==3)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','k','LineWidth',2);
        elseif(A(ii,3)==4)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','r','LineWidth',2);
        elseif(A(ii,3)==5)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','c','LineWidth',2);
        elseif(A(ii,3)==6)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','m','LineWidth',2);
        elseif(A(ii,3)==7)
        plot(center(1),center(2),'o','MarkerSize', 30,'MarkerEdgeColor','y','LineWidth',4);
        elseif(A(ii,3)==8)
        plot(center(1),center(2),'o','MarkerSize', 30,'MarkerEdgeColor','g','LineWidth',4);
        elseif(A(ii,3)==9)
        plot(center(1),center(2),'o','MarkerSize', 30,'MarkerEdgeColor','b','LineWidth',4);
        elseif(A(ii,3)==10)
        plot(center(1),center(2),'o','MarkerSize', 30,'MarkerEdgeColor','k','LineWidth',4);
        elseif(A(ii,3)==11)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','r','LineWidth',4);
        elseif(A(ii,3)==12)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','c','LineWidth',4);
        elseif(A(ii,3)==13)
        plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','m','LineWidth',4);
        end
% elseif(jj==1)'MarkerFaceColor',[.49 1 .63],
%         if(A(ii,3)==0)
%         plot(center(1),center(2),'+','MarkerSize', 20,'MarkerFaceColor','y','LineWidth',2);
%         elseif(A(ii,3)==1)
%         plot(center(1),center(2),'+','MarkerSize', 20,'MarkerFaceColor','g','LineWidth',2);
%         elseif(A(ii,3)==2)
%         plot(center(1),center(2),'+','MarkerSize', 20,'MarkerFaceColor','b','LineWidth',2);
%         elseif(A(ii,3)==3)
%         plot(center(1),center(2),'+','MarkerSize', 20,'MarkerFaceColor','k','LineWidth',2);
%         elseif(A(ii,3)==4)
%         plot(center(1),center(2),'+','MarkerSize', 20,'MarkerFaceColor','r','LineWidth',2);
%         end
%  elseif(jj==2)
%         if(A(ii,3)==0)
%         plot(center(1),center(2),'*','MarkerSize', 10,'MarkerFaceColor','y');
%         elseif(A(ii,3)==1)
%         plot(center(1),center(2),'*','MarkerSize', 10,'MarkerFaceColor','g');
%         elseif(A(ii,3)==2)
%         plot(center(1),center(2),'*','MarkerSize', 10,'MarkerFaceColor','b');
%         elseif(A(ii,3)==3)
%         plot(center(1),center(2),'*','MarkerSize', 10,'MarkerFaceColor','k');
%         elseif(A(ii,3)==4)
%         plot(center(1),center(2),'*','MarkerSize', 10,'MarkerFaceColor','r');
%         end
%  else
%         if(A(ii,3)==0)
%         plot(center(1),center(2),'>','MarkerSize', 10,'MarkerFaceColor','y');
%         elseif(A(ii,3)==1)
%         plot(center(1),center(2),'>','MarkerSize', 10,'MarkerFaceColor','g');
%         elseif(A(ii,3)==2)
%         plot(center(1),center(2),'>','MarkerSize', 10,'MarkerFaceColor','b');
%         elseif(A(ii,3)==3)
%         plot(center(1),center(2),'>','MarkerSize', 10,'MarkerFaceColor','k');
%         elseif(A(ii,3)==4)
%         plot(center(1),center(2),'>','MarkerSize', 10,'MarkerFaceColor','r');
%         end
%  end

end
end