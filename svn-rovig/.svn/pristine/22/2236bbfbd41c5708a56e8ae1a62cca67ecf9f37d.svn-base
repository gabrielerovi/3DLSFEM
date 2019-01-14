function maps=ArnoldMonotoneConstraint(mesh,maps,parameters)

L=length(mesh);

dim=parameters.dim;

toll=10^(-16);

for lev=1:L-1
    
    C=lev;
    F=lev+1;
    
    NC=mesh{C}.N;
    NEC=mesh{C}.NE
    NF=mesh{F}.N;
    nodeC=mesh{C}.node;
    nodeF=mesh{F}.node;
    edgeC=mesh{C}.edge;
    edgeF=mesh{F}.edge;
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
    
    
    
    for eC=1:NEC
        pC=edgeC(eC,:);
        pointC=0.5*(nodeC(pC(1),1:dim)+nodeC(pC(2),1:dim));
        
        for ii=1:2
        patchF_pC=maps.Patch_Edge{F}{pC(ii)};
        
        for eeF=patchF_pC
            pF=setdiff(edgeF(eeF,:),pC(ii));
            
            if(norm(nodeF(pF,1:dim)-pointC)<toll)
                maps.Patch_Edge_Monotone{C}{eC}(ii)=eeF;
            end
            
            
        end
        
        end
        
     

        
    end


end