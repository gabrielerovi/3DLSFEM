function mesh=mesh_contact_dofs(mesh,parameters)
L=length(mesh);
[dirichletN,n_and_or_tN,bool_bcN]= boundary_value_bool(1);
[dirichletE,n_and_or_tE,bool_bcE]= boundary_value_bool(2);


for lev=1:L
    
    N=mesh{lev}.N;
    NE=mesh{lev}.NE;
    N_bc=mesh{lev}.N_bc;
    E_bc=mesh{lev}.E_bc;
    
    contN=0;
    contE=0;
    mesh{lev}.N_bc_contact=sparse(N,1);
    mesh{lev}.E_bc_contact=sparse(NE,1);
    mesh{lev}.N_contact=[];
    mesh{lev}.E_contact=[];
    for nn=1:N
        if(N_bc(nn)>0 &&dirichletN(N_bc(nn),3)==1 )
            contN=contN+1;
            mesh{lev}.N_bc_contact(nn)=nn;
            mesh{lev}.N_contact(contN)=nn;
        end
    end
    
    for ee=1:NE
         if(E_bc(ee)>0 &&dirichletE(E_bc(ee),3)==1 )
            contE=contE+1;
            mesh{lev}.E_bc_contact(ee)=ee;
            mesh{lev}.E_contact(contE)=ee;
        end
    end  
    
%     if(~isempty(mesh{lev}.N_contact))
%     mesh{lev}.N_contact_map=containers.Map(mesh{lev}.N_contact,1:length(mesh{lev}.N_contact));
%     end
%     if(~isempty(mesh{lev}.E_contact))
%     mesh{lev}.E_contact_map=containers.Map(mesh{lev}.E_contact,1:length(mesh{lev}.E_contact));
%     end
     
   
 mesh{lev}.Remove=[mesh{lev}.E_remove, mesh{lev}.E_remove + mesh{lev}.NE, mesh{lev}.N_remove + 2 * mesh{lev}.NE, mesh{lev}.N_remove + 2 * mesh{lev}.NE + mesh{lev}.N ];
 mesh{lev}.RemoveNT=unique([mesh{lev}.Remove,mesh{lev}.E_contact + mesh{lev}.NE ]);
 mesh{lev} = rmfield(mesh{lev},'Remove');
end


end