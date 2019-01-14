


function [AArnoldLocal]=ArnoldLocalAssembling(lev,mesh,maps,parameters)

householder=1;
[M_Normal_Tangent,M_Normal_TangentT] = MatrixOnGammaCwithNormalTangentComponents(mesh{lev},householder);


N=mesh{lev}.N;
NE=mesh{lev}.NE;
N_to_T=mesh{lev}.N_to_T;
qrule=parameters.qrule;
alpha=parameters.alpha;
beta=parameters.beta;

N_dirichlet=mesh{lev}.N_dirichlet;
E_dirichlet=mesh{lev}.E_dirichlet;
E_contact=mesh{lev}.E_contact;



C_eq=parameters.C_eq;
C_const=parameters.C_const;
C_asym=parameters.C_asym;
input_name=parameters.input_name;
meshLoc.edge_per_elem=mesh{lev}.edge_per_elem;
meshLoc.node_per_elem=mesh{lev}.node_per_elem;

EmapGlob2Loc=maps.EmapGlob2Loc{lev};
NmapGlob2Loc=maps.NmapGlob2Loc{lev};
Patch_Node=maps.Patch_Node{lev};
Patch_Edge=maps.Patch_Edge{lev};

for nn=1:N
    
    N_Glob=Patch_Node{nn};
    E_Glob=Patch_Edge{nn};
    
    E_Glob1=E_Glob;
    E_Glob2=E_Glob+NE;
    N_Glob1=N_Glob+2*NE;
    N_Glob2=N_Glob1+N;  
    Tot_Glob=[E_Glob1,E_Glob2,N_Glob1,N_Glob2];
    M_NT_Loc=M_Normal_Tangent(Tot_Glob,Tot_Glob);
    
    nodes=maps.Patch_Node{lev}{nn};
    nodes_bc=maps.Patch_Boundary_Node{lev}{nn};
    nodes_bc=cell2mat(values(NmapGlob2Loc{nn},num2cell(nodes_bc,1)));  
    edges=maps.Patch_Edge{lev}{nn};
    elems=cell2mat(N_to_T{nn});
    elemE=mesh{lev}.elemE(elems,:);
    
    
    meshLoc.node=mesh{lev}.node(nodes,:);
    meshLoc.edge=mesh{lev}.edge(edges,:);
    meshLoc.elem=mesh{lev}.elem(elems,:);
    
    
    for tt=1:length(elems)
    meshLoc.elem(tt,:)=cell2mat(values(NmapGlob2Loc{nn},num2cell(meshLoc.elem(tt,:),1)));  
    meshLoc.elemE(tt,:)=cell2mat(values(EmapGlob2Loc{nn},num2cell(elemE(tt,:),1)));  
    end
    
    meshLoc.edge=cell2mat(values(EmapGlob2Loc{nn},num2cell(edges,1)));
    
    meshLoc.N=length(nodes);
    meshLoc.NE=length(edges);
    meshLoc.NT=length(elems);
    
    % (1,1),(1,2),(2,2)
    for ii=1:2
        for jj=ii:2
             ALoc{ii,jj}=assembling2DRTRT(meshLoc,qrule,alpha,beta,ii,jj,C_eq,C_const,C_asym,input_name);
        end
    end

    % (1,3),(1,4), (2,3),(2,4)
    for ii=1:2
        for jj=3:4
             ALoc{ii,jj}=assembling2DRTP1(meshLoc,qrule,alpha,beta,ii,jj,C_eq,C_const,input_name);
        end
    end
    % (3,1),(3,2), (4,1),(4,2)
    for ii=3:4
        for jj=1:2
             ALoc{ii,jj}=assembling2DRTP1(meshLoc,qrule,alpha,beta,ii,jj,C_eq,C_const,input_name);
        end
    end
    
    % (3,3),(3,4),(4,4)
    for ii=3:4
        for jj=ii:4
             ALoc{ii,jj}=assembling2DP1P1(meshLoc,qrule,alpha,beta,ii,jj,C_eq,C_const,input_name);
        end
    end
    
    
    AArnoldLocal{nn}=[ALoc{1,1}  ALoc{1,2} ALoc{1,3}  ALoc{1,4};
                      ALoc{1,2}' ALoc{2,2} ALoc{2,3}  ALoc{2,4};
                      ALoc{3,1}  ALoc{3,2} ALoc{3,3}  ALoc{3,4};
                      ALoc{4,1}  ALoc{4,2} ALoc{3,4}' ALoc{4,4};        ];
    
                  
   % go into the nt system               
   AArnoldLocal{nn} = M_NT_Loc *AArnoldLocal{nn} * M_NT_Loc';
   
   % remove boundary nodes, sicne they are bcs
   bc=2*meshLoc.NE+nodes_bc;
   bc=[bc,bc+meshLoc.N];
   AArnoldLocal{nn}(bc,:)=0;
   for ii=bc
       AArnoldLocal{nn}(ii,ii)=1;
   end
   
   % remove real bcs
    N_bc=find(N_dirichlet(N_Glob)==1);
    E_bc=find(E_dirichlet(E_Glob)==1);
    remove=[E_bc,E_bc+meshLoc.NE, N_bc+2*meshLoc.NE,N_bc+2*meshLoc.NE+meshLoc.N];
    
    
    AArnoldLocal{nn}(remove,:)=0;
    for rr=remove
        AArnoldLocal{nn}(rr,rr)=1;
    end
   % remove frictionless edges
   E_contact_Loc=intersect(E_Glob,E_contact);
   if(~isempty(E_contact_Loc))
   E_contact_Loc=cell2mat(values(EmapGlob2Loc{nn},num2cell(E_contact_Loc,1)));
   
   E_contact_Loc=E_contact_Loc+meshLoc.NE;
   AArnoldLocal{nn}(E_contact_Loc,:)=0;
    for rr=E_contact_Loc
        AArnoldLocal{nn}(rr,rr)=1;
    end   
   end
end

end

