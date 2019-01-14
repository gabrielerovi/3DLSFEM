function mesh=mesh_contact_dofs(mesh,parameters)
L=length(mesh);
[dirichletN,n_and_or_tN,bool_bcN]= boundary_value_bool(1);
[dirichletF,n_and_or_tF,bool_bcF]= boundary_value_bool(3);


for lev=1:L
    
    N=mesh{lev}.N;
    NF=mesh{lev}.NF;
    N_bc=mesh{lev}.N_bc;
    F_bc=mesh{lev}.F_bc;
    
    contN=0;
    contF=0;
    mesh{lev}.N_bc_contact=sparse(N,1);
    mesh{lev}.F_bc_contact=sparse(NF,1);
    mesh{lev}.N_contact=[];
    mesh{lev}.F_contact=[];
    for nn=1:N
        if(N_bc(nn)>0 &&dirichletN(N_bc(nn),end)==1 )
            contN=contN+1;
            mesh{lev}.N_bc_contact(nn)=nn;
            mesh{lev}.N_contact(contN)=nn;
        end
    end
    
    for ff=1:NF

         if(F_bc(ff)>0 &&dirichletF(F_bc(ff),end)==1 )
            contF=contF+1;
            mesh{lev}.F_bc_contact(ff)=ff;
            mesh{lev}.F_contact(contF)=ff;
        end
    end
    
    mesh{lev}.Remove=[mesh{lev}.F_remove, ...
                   mesh{lev}.F_remove + mesh{lev}.NF, ...
                   mesh{lev}.F_remove + 2 * mesh{lev}.NF, ...
                   mesh{lev}.N_remove + 3 * mesh{lev}.NF, ...
                   mesh{lev}.N_remove + 3 * mesh{lev}.NF + mesh{lev}.N, ... 
                   mesh{lev}.N_remove + 3 * mesh{lev}.NF + 2* mesh{lev}.N ];
 mesh{lev}.RemoveNT=unique([mesh{lev}.Remove,mesh{lev}.F_contact + mesh{lev}.NF , mesh{lev}.F_contact + mesh{lev}.NF * 2 ]);
end


end