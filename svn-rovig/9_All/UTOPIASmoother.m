
function [x]=UTOPIASmoother(A,x,b,vertices_per_processor,Patch_Internal_All,smoothing_steps) 


num_procs=length(vertices_per_processor);
%  SMOOTHING-STEPS

res=b- A * x;
% REMARK: this methods only works if each processor has both nodal and face dofs of the same patch
% LIMIT CASE: processor 1 has only nodes, processor 2 has only face dofs -> face-dofs never updated!

% Now we solve A correction = res, on each processor
% We loop on all the processor pp
% Then on all the nodes of the given processor pp we do a patch solve
% We locally update only the values belonging to the given processor
% Other values are set to zero as bcs

% correction initially set to zero
correction=zeros(length(x),1);
for jj=1:smoothing_steps
for pp=1:num_procs
for nn=vertices_per_processor{pp}

    % define dofs patch
    Patch=Patch_Internal_All{nn};
    % build the local matrix
    A_Loc= A (Patch,Patch);
    % find local bcs, i.e. dofs not belonging to the current processor
    bc_dofs_Loc=find(dofs_inside_processor{pp}(Patch)==0);
    % set local bcs on A_loc and b_loc
    A_Loc(bc_dofs_Loc,:)=0;
    for ii=bc_dofs_Loc
    	A_Loc(ii,ii)=1;
        % we are not going to modify values on other processors, thus correction=0 on these dofs
        b_Loc(ii)=0.0;
    end

    b_Loc= res (Patch)- A(Patch,:) * x;
    c_Loc=A_Loc\b_Loc;

    x(Patch)=x(Patch) + c_Loc;
    
    b_Loc= b (Patch)- A(Patch,:) * correction;
    c_Loc=A_Loc\b_Loc;
    correction(Patch)=correction(Patch) + c_Loc;
end  
end
end


% synchronize x, knowing that dofs in correction are not shared among processors,
% therefore: no need for average or other tricks
x=x+correction;


end
























