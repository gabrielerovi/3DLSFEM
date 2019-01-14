function n_cont=normal_contact(node,parameters)

if(parameters.dim==2)
     n_cont=[ parameters.body_normal_x( node(1), node(2));
              parameters.body_normal_y( node(1),  node(2));];
else
     n_cont=[ parameters.body_normal_x( node(1), node(2),node(3));
              parameters.body_normal_y( node(1), node(2),node(3));
              parameters.body_normal_z( node(1), node(2),node(3));];   
end
end