function [x,is_constrained]=ContactConstraint(mesh,x,dof,type_of_dof,parameters)

% x is a 2d vector [ux,uy] or [(sigma n) x, (sigma n) y]
% dof is the node-edge position for the first component.
% if type_of_dof==1 (node), then 1<= dof <= N
% if type_of_dof==2 (edge), then 1<= dof <= NE
% from dof I know the node-edge position in mesh


is_constrained=0;

if(type_of_dof==1)
    % position
    position=mesh.node(dof,:);
    % gap function
    gap=parameters.gap;
    % HH transformation to describe x in local coordinates:
    % the new first component is directed as the normal
    normal=mesh.normal_node{dof};
%     normal_patch=mesh.normal_node_patch{dof};
%     dotprod1=normal_patch(:,1)'*x;
%     dotprod2=normal_patch(:,2)'*x;
%     if(dotprod1>dotprod2)
%         normal=normal_patch(:,1);
%     else
%         normal=normal_patch(:,2);
%     end
    [xn,H]=HouseHolderTransformation(x,normal);

    if(xn(1)>=gap(position(1),position(2)))
     is_constrained=1;
     xn(1)=gap(position(1),position(2));
    end
    x=H'*xn;
else 
        % tangent of the edge
        normal_edge=mesh.normal_edge{dof};
        % normal of the edge
        tangent_edge=[normal_edge(2); -normal_edge(1)];
        
        % phidotn= phi cdot n
        phidotn=phi_dot_n(mesh,dof);
        % (sigma n)= (phi n) Coeff
        normalstress=x*phidotn;
        
        % in case of frictionless contact, we remove the tangent component
        if(parameters.frictionless==1)
        [normalstress_tangentdirection,H]=HouseHolderTransformation(normalstress,tangent_edge);
        normalstress_tangentdirection(1)=0;
        normalstress=H'*normalstress_tangentdirection; 
        end
        
        % then, we enforce the negative pressure condition
        [normalstress_normaldirection,H]=HouseHolderTransformation(normalstress,normal_edge);
        if(normalstress_normaldirection(1)>=0)
            is_constrained=1;
            normalstress_normaldirection(1)=0;
        end
    % then we project it back to the standard basis    
    normalstress=H'*normalstress_normaldirection;      
    % now we want to recompute the coefficients [phi n C1, phi n C2]=[x1,x2],
    % i.e. C1=x1/ (phi *n), C2=x2/ (phi *n)
    % remark: the shape function is the same, so we just have to compute
    % only phi * n
    % coefficients C1, C2
    x=normalstress./phidotn;
    
end


end