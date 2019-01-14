function maps=ArnoldMonotoneConstraint(mesh,maps,parameters)

L=length(mesh);

dim=parameters.dim;

toll=10^(-16);
maps.Patch_Face_Monotone=[];
maps.Patch_Node_Monotone=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop on all the coarse level C and for each level:                        %%%%%%%%%%%%%%
%%% Loop on al the coarse nodes nC. Find the coarse and fine patches PC,PF    %%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for lev=1:L-1
    
    C=lev;
    F=lev+1;
    
    NC=mesh{C}.N;
    NFC=mesh{C}.NF;
    NF=mesh{F}.N;
    nodeC=mesh{C}.node;
    nodeF=mesh{F}.node;
    faceC=mesh{C}.face;
    faceF=mesh{F}.face;
    % loop on all coarse nodes
    for nC=1:NC
        
        % define the coarse patch nodes, excluding the actual one nC
        patchC=setdiff(maps.Patch_Node{C}{nC},nC);
        
        patchF_nC=setdiff(maps.Patch_Node{F}{nC},nC);
        
        % loop on all the coarse patch nodes
        cont_patchC=0;
        for pnC=patchC
            
            cont_patchC=cont_patchC+1;
            
            % for each coarse patch node, search the fine patch
            patchF_pnC=setdiff(maps.Patch_Node{F}{pnC},pnC);
            
            patchF_intersect=intersect(patchF_nC,patchF_pnC);
            
            pointC=0.5 * (nodeC(nC,1:dim)+nodeC(pnC,1:dim));
            
            for ii =patchF_intersect
                
                pointF=nodeF(ii,1:dim);
                if(norm(pointC-pointF)<toll)
                  maps.Patch_Node_Monotone{C}{nC}{cont_patchC}=[pnC;ii];
                break
                end
            end
            
        end
        
        
    end
    


end


for lev=1:L-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Loop on all coarse faces                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C=lev;
    F=lev+1;
    
    NC=mesh{C}.N;
    NFC=mesh{C}.NF;
    NF=mesh{F}.N;
    nodeC=mesh{C}.node;
    nodeF=mesh{F}.node;
    faceC=mesh{C}.face;
    faceF=mesh{F}.face;

       
    for fC=1:NFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Find the midpoint nodes of the face F (node n1,n2,n3) by:      %%%%%%%
%%%% For each n = n1,n2,n3:                                         %%%%%%%
%%%% Loop on the Patch_Node_Monotone to find edges=[dofC,dofF];     %%%%%%%
%%%% Then check if dofF is one of the midpoint nodes                %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pC=faceC(fC,:);
        nodepC=nodeC(pC,1:dim);
        pointC(1:dim,1)=0.5 * (nodepC(1,:)+nodepC(2,:) )';
        pointC(1:dim,2)=0.5 * (nodepC(2,:)+nodepC(3,:) )';
        pointC(1:dim,3)=0.5 * (nodepC(1,:)+nodepC(3,:) )';        
        RelMatrix=[0 1 3;1 0 2;3 2 0];
        midpointnodes=zeros(3,1);
        for ii=1:dim
            for kk1=1:length(maps.Patch_Node_Monotone{C}{pC(ii)})
                dofF=maps.Patch_Node_Monotone{C}{pC(ii)}{kk1}(2);
                point_patchF=nodeF(dofF,1:dim)';  
                
                for kk2=1:dim
                    if(norm(pointC(:,kk2)-point_patchF)<toll && RelMatrix(ii,kk2)>0)
                        midpointnodes(RelMatrix(ii,kk2),1)=dofF;
                    end
                end
                
            end
            
        end    
            
   
    
    coarsefacesF=[pC(1)           ,midpointnodes(1), midpointnodes(3) ;
                  pC(2)           ,midpointnodes(1), midpointnodes(2) ;
                  pC(3)           ,midpointnodes(2), midpointnodes(3) ;
                  midpointnodes(1),midpointnodes(2), midpointnodes(3);];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Once the midpoint nodes of the face F are known, search the fine faces %%%%%%
%%%% Loop on 2 of the midpoint nodes using Patch_Face_Monotone(F)           %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
         for ii=1:dim-1
             pF=midpointnodes(ii);
             fine_faces_patch=maps.Patch_Face{F}{pF};
             
             for kk1=1:length(fine_faces_patch)
                 
                 for kk2=1:dim+1
                     if(isequal(sort(coarsefacesF(kk2,:)),sort(faceF(fine_faces_patch(kk1),:))))
                         maps.Patch_Face_Monotone{C}{fC}(kk2)=fine_faces_patch(kk1);
                     end                                          
                 end
             end
             
             
             
         end
%         patchF_pC=maps.Patch_Face{F}{pC(ii)};
%         
%         for ffF=patchF_pC
%             
%             barycenterPatchF_pC=mean(nodeF(faceF(ffF,:),:));
%             
%             if(norm(nodeF(barycenterPatchF_pC,1:dim)-barycenter(kk))<toll)
%                 maps.Patch_Face_Monotone{C}{fC}(ii)=ffF;
%             end
%             
%             
%             pF=setdiff(faceF(ffF,:),pC(ii));
%             
%             if(length(pF)<length(faceF(ffF,:)))
%             if(norm(nodeF(pF,1:dim)-pointC)<toll)
%                 maps.Patch_Face_Monotone{C}{fC}(ii)=ffF;
%             end
%             end
%             
%         end
%         
%         end
        
     

        
    end 
 
 end

end
    
    
    