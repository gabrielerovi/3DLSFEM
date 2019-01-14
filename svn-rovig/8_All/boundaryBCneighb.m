function [NeumannNode,BCneighb,BCneighbCoeff]=boundaryBCneighb(mesh)

E_remove=mesh.E_remove;
NeumannEdge=mesh.edge(E_remove,:);
NeumannEdgeLength=length(NeumannEdge(:,1));
temp=unique(mesh.edge(mesh.E_remove,:));
N_dirichlet=mesh.N_dirichlet;
E_to_T=mesh.E_to_T(mesh.E_remove);

cont=0;
NeumannNode=[];
NeumannNodeMap=zeros(mesh.N,1);
for i=1:length(temp)
    if(N_dirichlet(temp(i))==0)
       cont=cont+1;
       NeumannNode(cont)=temp(i);
       NeumannNodeMap(temp(i))=cont;
    end
end
if(isempty(NeumannNode))
    NeumannNodeLength=0;
else
NeumannNodeLength=length(NeumannNode);
end
BCneighb=zeros(NeumannNodeLength,2);
BCneighbCoeff=cell(NeumannNodeLength,1);

for ee=1:NeumannEdgeLength
    
    edge=E_remove(ee);
    T=cell2mat(E_to_T{ee});
    elem=mesh.elem(T,:);
    
    nodes=NeumannEdge(ee,:);
    node=mesh.node(nodes,:);
    extranode=setdiff(elem, nodes);
    nodes=[nodes extranode];
    node=[node;mesh.node(extranode,:) ];
    grad = P1grad2D  ( node  );
    
    
    
    temp=(0.5*(node(1,:)+node(2,:))-node(3,:) );
    normal=[node(1,2)-node(2,2) ; -(node(1,1)-node(2,1))];
    sign=temp*normal;
    sign=sign/abs(sign); 
    
    normal= sign*normal/norm(normal)  ;
    coeff1=grad(1,2)*normal(1) - grad(1,1)* normal(2);
    coeff2=grad(2,2)*normal(1) - grad(2,1)* normal(2);
    
    NN1=NeumannNodeMap(nodes(1));
    NN2=NeumannNodeMap(nodes(2));
    if(N_dirichlet(nodes(1))==0 && NN1>0)
        if(isempty(BCneighbCoeff{NN1}))
       BCneighbCoeff{NN1}=coeff1;
       BCneighbCoeff{NN1}=[BCneighbCoeff{NN1},coeff2];
        else
       BCneighbCoeff{NN1}=BCneighbCoeff{NN1}+coeff1;
       BCneighbCoeff{NN1}=[BCneighbCoeff{NN1},coeff2];
        end
        if(BCneighb(NN1,1)==0)
       BCneighb(NN1,1)=nodes(2); 
        else
       BCneighb(NN1,2)=nodes(2); 
        end
    end
     if(N_dirichlet(nodes(2))==0 && NN2>0)
        if(isempty(BCneighbCoeff{NN2}))
       BCneighbCoeff{NN2}=coeff2;
       BCneighbCoeff{NN2}=[BCneighbCoeff{NN2},coeff1];
        else
       BCneighbCoeff{NN2}(1)=BCneighbCoeff{NN2}(1)+coeff2;
       BCneighbCoeff{NN2}=[BCneighbCoeff{NN2},coeff1];
        end       
       
        if(BCneighb(NN2,1)==0)
       BCneighb(NN2,1)=nodes(1); 
        else
       BCneighb(NN2,2)=nodes(1); 
        end 
    end   

        
    end
    
    
    
    
       
end

