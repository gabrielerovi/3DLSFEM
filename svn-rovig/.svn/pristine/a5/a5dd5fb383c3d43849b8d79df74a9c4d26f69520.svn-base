function M=assembling_UU(qrule,node,node_per_elem,alpha,beta,Ceq,Cconst,Casym)

Volume=VolumeTetrahedronAndNormalsigns(node) ;
[P1grad] = P1grad3D(node);
M=cell(3,3);
for mm=1:3
for nn=mm:3
M{mm,nn}=zeros(node_per_elem,node_per_elem);
end
end


   for n1=1:node_per_elem    
        for n2=1:node_per_elem
                        

            value(1,1)= Volume * Cconst * (      P1grad(n1,1)*P1grad(n2,1) + 0.5 * P1grad(n1,2)*P1grad(n2,2) + 0.5 * P1grad(n1,3)*P1grad(n2,3));  
                 
            value(1,2)= Volume * Cconst * (0.5 * P1grad(n1,2)*P1grad(n2,1));

            value(1,3)= Volume * Cconst * (0.5 * P1grad(n1,3)*P1grad(n2,1));
                 
            value(2,2)= Volume * Cconst * (0.5 * P1grad(n1,1)*P1grad(n2,1) + 	   P1grad(n1,2)*P1grad(n2,2) + 0.5 * P1grad(n1,3)*P1grad(n2,3));  
                 
            value(2,3)= Volume * Cconst * (0.5 * P1grad(n1,3)*P1grad(n2,2));

            value(3,3)= Volume * Cconst * (0.5 * P1grad(n1,1)*P1grad(n2,1) + 0.5 * P1grad(n1,2)*P1grad(n2,2) +        P1grad(n1,3)*P1grad(n2,3));  
                 
            
          
          for mm=1:3
              for nn=mm:3
                   M{mm,nn}(n1,n2)=M{mm,nn}(n1,n2)+value(mm,nn);  
              end
          end
          
        end
    end

          for mm=1:2
              for nn=mm+1:3
                   M{nn,mm}=M{mm,nn}';
              end
          end


end
