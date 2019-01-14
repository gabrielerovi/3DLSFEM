function [edge_boundary_new,reorder,lengths]=nodes_on_edge_boundary(mesh)

% vector of label of the E_bc
E_remove=mesh.E_remove;
% edges
edge=mesh.edge(E_remove,:); 
if(~isempty(edge))
% vector of label of the E_bc
bc_label=mesh.E_bc(mesh.E_remove);


% number of dirichlet edge
NER=length(E_remove);

% now in bc label we have the number of the boundaries. we want to map this
% into 1...numberoflabels
u=unique(bc_label);
numberoflabel=length(u);


for nn1=1:NER

    for nn2=1:numberoflabel
        if(bc_label(nn1)==u(nn2))
            bc_label(nn1)=nn2;
        end
    end
        
end




% number of -edge dirichlet edge
Nbc=length(unique(bc_label));
edge_boundary=cell(Nbc,1);
for nn=1:Nbc
    edge_boundary{nn}=[];
end

for ee=1:NER
     
    nn=bc_label(ee);
    
    edge_boundary{nn}(end+1)=edge(ee,1);
    edge_boundary{nn}(end+1)=edge(ee,2);
    
end

for nn=1:Nbc
    edge_boundary{nn}=unique(sort(edge_boundary{nn}));
   % first(nn)=edge_boundary{nn}(1);
end

intersection=zeros(Nbc,Nbc);
edge_boundary_new=cell(1,1);
real_bc=[];
cont=0;
if(Nbc>1)
 cont=0;   
for nn1=1:Nbc
    for nn2=nn1+1:Nbc
        if( ~isempty(intersect(edge_boundary{nn1},edge_boundary{nn2}) ) )
            intersection(nn1,nn2)=1;
        end
    end
    
    
    
end
remove=[];
for nn1=1:Nbc
    edge_boundary_new{nn1}=edge_boundary{nn1};
    for nn2=nn1+1:Nbc
        if( intersection(nn1,nn2)) 
            cont=cont+1;
            remove(cont)=nn2;
            edge_boundary_new{nn1}=cat(2,edge_boundary_new{nn1},edge_boundary{nn2});
            for nn3=nn1+1:Nbc
            if( intersection(nn2,nn3))
            cont=cont+1;
            remove(cont)=nn3;
            edge_boundary_new{nn1}=cat(2,edge_boundary_new{nn1},edge_boundary{nn3});
            end
            
            end
        end
        
    end   
end




remove=unique(remove);
edge_boundary_new(remove) = []   
for ii=1:length(edge_boundary_new)
    edge_boundary_new{ii}=unique(edge_boundary_new{ii});
end
else
    edge_boundary_new{1}=edge_boundary{1};
end



vv=1:mesh.N;
for ii=1:length(edge_boundary_new)
   vv=setdiff(vv,edge_boundary_new{ii});
end

reorder=[];
for ii=1:length(edge_boundary_new)
   reorder=cat(2,reorder,edge_boundary_new{ii});
   lengths(ii)=length(edge_boundary_new{ii});
end
reorder=cat(2,reorder,vv);

lengths(ii+1)=length(vv);

else
  edge_boundary_new=[];
  reorder=1:mesh.N;
  lengths(1)=0;
  lengths(2)=mesh.N;
end




end


