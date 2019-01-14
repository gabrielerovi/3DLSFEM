function duality_complementarity=ContactIntegral(mesh,x,parameters)

gap=parameters.gap;
L=size(mesh);
L=L(1);

E_bc=mesh{L}.E_bc;
NE=mesh{L}.NE;
N=mesh{L}.N;
type_of_dof=2;
[dirichlet,n_and_or_t,bool_bc]= boundary_value_bool(type_of_dof);

sigma=x(1:2*NE);
disp=x(1+2*NE:end);

duality_complementarity=0;

for ee=1:NE
    label=E_bc(ee);
    % on boundary
    if(label>0)
        % in contact
        if(dirichlet(label,3)==1)
        edge=mesh{L}.edge(ee,:);
        node1=mesh{L}.node(edge(1),[1,2]);
        node2=mesh{L}.node(edge(2),[1,2]);
        
        disp_loc1=disp([edge(1),edge(1)+N]);
        disp_loc2=disp([edge(2),edge(2)+N]);

        [disp_loc1,Nconstraint1]=ContactConstraint(mesh{L},disp_loc1,edge(1),1,parameters);
        [disp_loc2,Nconstraint2]=ContactConstraint(mesh{L},disp_loc2,edge(2),1,parameters);

        
        gap1=gap(node1(1),node1(2));
        gap2=gap(node2(1),node2(2));
        
        phidotn=phi_dot_n(mesh{L},ee);
        normalstress=sigma([ee,ee+NE])*phidotn;
        [normalstress_normaldirection,H]=HouseHolderTransformation(normalstress,mesh{L}.normal_edge{ee});            
            
        
        duality_complementarity= duality_complementarity + 0.5 * norm(node1-node2) * (( disp_loc1(1)-gap1)+ (disp_loc2(1)-gap2)) * normalstress_normaldirection(1) ;
       % [ee,normalstress_normaldirection(1), ( disp_loc1(1)-gap1), (disp_loc2(1)-gap2),duality_complementarity]
        end
    end
end

end