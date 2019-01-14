
function print_stress_error(mesh,x,component,exact)

L=size(mesh);
L=L(1);
figure

set(gca,'FontSize',20)
txt1='FLUX ';
if(component==1)
    txt2=' X';
else
    txt2=' Y';
end

txt=strcat(txt1,txt2);
title(txt);
h=colorbar('northoutside');
colormap('jet');
hold on
for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    elemE=mesh{L}.elemE(t,:);
    node=mesh{L}.node(elem,:);
    line(node([1:end,1],1),node([1:end,1],2));
end
minimum=0;
maximum=0;
for t=1:mesh{L}.NT
    elem=mesh{L}.elem(t,:);
    elemE=mesh{L}.elemE(t,:);
    node=mesh{L}.node(elem,:);
    v=x(elemE);
    phi=phiRT2D(node,node);
    vector=zeros(length(phi(:,1)),2);
%     nodeE= [ 0.5 * ( node(1,1)+node(2,1) ),0.5 * (node(1,2)+node(2,2));
%              0.5 * ( node(2,1)+node(3,1) ),0.5 * (node(2,2)+node(3,2));
%              0.5 * ( node(1,1)+node(3,1) ),0.5 * (node(1,2)+node(3,2));];
    
    barycenter=sum(node)/length(node(:,1));
    for e=1:length(phi(:,1))
        for i=1:3
            for dd=1:2
            vector(e,dd)= vector(e,dd)+v(i)*phi(i,e,dd);
            end
        end
    end

    for ii=1:length(node(:,1))
     exe(ii,1)=exact(node(ii,1),node(ii,2));
    end
    
     vector=abs(exe-vector);
     www=[vector(1,component);vector(2,component);vector(3,component)];
     
     mean=1./3.0 *(www(1)+www(2)+www(3));
     www=[mean;mean;mean];
     fill3(node(:,1),node(:,2),www,www);
     
    minimum=(min(minimum,min(www)));
    maximum=(max(maximum,max(www)));
end


set(gca,'FontSize',20)
txt1='ERROR STRESS ';

set(h,'YTick',[minimum:1:maximum]);

end





