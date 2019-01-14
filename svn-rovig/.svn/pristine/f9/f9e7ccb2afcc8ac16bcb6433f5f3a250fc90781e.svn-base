function p = ContactPressure(mesh,xnt,parameters)




L=size(mesh);
L=L(1);

N=mesh{L}.N;
NF=mesh{L}.NF;


node=mesh{L}.node;
dim=length(node(1,:));
face=mesh{L}.face;

% the first component is the pressure
sigman=xnt(1:NF);


%initialize contact pressure
nodecontact=sparse(N,1);
pressurenodecontact=sparse(N,1);
for ff=mesh{L}.F_contact  
       nodeface=node(face(ff,1:dim),1:dim);
       nodecontact(face(ff,1:dim))=nodecontact(face(ff,1:dim))+1;
       pressurenodecontact(face(ff,1:dim))=pressurenodecontact(face(ff,1:dim))+sigman(ff);
end

for nn=1:mesh{L}.N
    if(nodecontact(nn)>0)
    pressurenodecontact(nn)=pressurenodecontact(nn)/nodecontact(nn);    
    end
end
figure
hold on

% for ff=mesh{L}.F_contact  
%        nodeface=node(face(ff,1:dim),1:dim);
%        fill3(nodeface(:,1),nodeface(:,2),pressurenodecontact(face(ff,1:dim)),pressurenodecontact(face(ff,1:dim)))
% end

l.FontSize=46;
xx=xlabel('X');
yy=ylabel('Y');
zz=zlabel('P');
set(gca,'fontsize',48);
title('PRESSURE');

figure
hold on
for ff=mesh{L}.F_contact  
    
       nodeface=node(face(ff,1:dim),1:dim);
       height=ones(3,1)*sigman(ff);
       fill3(nodeface(:,1),nodeface(:,2),height,height)

end

end