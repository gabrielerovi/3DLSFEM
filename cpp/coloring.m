close all
%clear all
%close all
s1='output/coloring_';
s3='.txt';
sdof1='output/dof';
num_procs=4;

for jj=0:num_procs-1
s2=num2str(jj);
s=strcat(s1,s2,s3);
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
xlim([0,1]);
ylim([0,1]);
for ii=1:n
center=[A(ii,1),A(ii,2)];
% if(jj==0)
txt = strcat(num2str(A(ii,5)),', ',num2str(A(ii,4)),',use=',num2str(A(ii,6)));
if(A(ii,5)<A(ii,7))
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
end
end
end





for jj=0:num_procs-1
s2=num2str(jj);
s=strcat(s1,s2,s3);
A=load(s);   
n=length(A(:,1));
 figure
hold on

for ii=1:n
center=[A(ii,1),A(ii,2)];
% if(jj==0)
txt = num2str(A(ii,6));
text(center(1),center(2),txt)
if(A(ii,6)==1)
plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','y','LineWidth',2);
end

end
end


 figure
hold on
for jj=0:num_procs-1
s2=num2str(jj);
s=strcat(s1,s2,s3);
A=load(s);   
n=length(A(:,1));


for ii=1:n

if(A(ii,6)==1)
plot(center(1),center(2),'o','MarkerSize', 20,'MarkerEdgeColor','y','LineWidth',2);
center=[A(ii,1),A(ii,2)];
% if(jj==0)
txt = num2str(A(ii,6));
text(center(1),center(2),txt)
end

end
end

