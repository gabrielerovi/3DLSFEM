function vertices_per_processor= UTOPIA_processsors(mesh,parameters,num_procs)

L=length(mesh);
 for lev=1:L
% compute the minimum and maximum x, y, then consider the bounding box and subdivide the mesh wrt this
min_coord=min(mesh{lev}.node);
max_coord=max(mesh{lev}.node);



sqr_subdivision=floor(sqrt(num_procs));
LX=max_coord(1)-min_coord(1);
LY=max_coord(2)-min_coord(2);
%then we make the residual subdivision in X
if(LX>LY)
DeltaX=LX/(1+sqr_subdivision);
DeltaY=LY/sqr_subdivision;
else
%then we make the residual subdivision in Y
DeltaX=LX/sqr_subdivision;
DeltaY=LY/(1+sqr_subdivision);
end


for pp1=1:sqr_subdivision
    for pp2=1:sqr_subdivision-1
        actual_procs=(pp1-1)*sqr_subdivision+pp2;                  
        bounding_box{actual_procs}(2,1)=(pp2-1)*DeltaX;
        bounding_box{actual_procs}(2,2)=(pp1-1)*DeltaY;
        bounding_box{actual_procs}(1,1)=pp2*DeltaX;
        bounding_box{actual_procs}(1,2)=pp1*DeltaY;
    end
       
        pp2=sqr_subdivision;
        actual_procs=(pp1-1)*sqr_subdivision+pp2;                  
        bounding_box{actual_procs}(2,1)=(pp2-1)*DeltaX+parameters.toll;
        bounding_box{actual_procs}(2,2)=(pp1-1)*DeltaY;
        bounding_box{actual_procs}(1,1)=pp2*DeltaX+parameters.toll;
        bounding_box{actual_procs}(1,2)=pp1*DeltaY;
end




sqr_subdivision_square=sqr_subdivision^2;
for pp=1+sqr_subdivision_square:num_procs
    if(LX>LY)
        pp2=1+sqr_subdivision;
        pp1=pp-sqr_subdivision_square;
    else
        pp2=pp-sqr_subdivision_square;
        pp1=1+sqr_subdivision;
    end
    if(pp==num_procs)
        bounding_box{pp}(2,1)=(pp2-1)*DeltaX+parameters.toll;
        bounding_box{pp}(2,2)=(pp1-1)*DeltaY;
        bounding_box{pp}(1,1)=pp2*DeltaX+parameters.toll;
        bounding_box{pp}(1,2)=pp1*DeltaY;  
    else
        bounding_box{pp}(2,1)=(pp2-1)*DeltaX;
        bounding_box{pp}(2,2)=(pp1-1)*DeltaY;
        bounding_box{pp}(1,1)=pp2*DeltaX;
        bounding_box{pp}(1,2)=pp1*DeltaY;
    end   
  end


% subdivide the dofs per processor
     for pp=1:num_procs
         vertices_per_processor{lev}{pp}=[];
%          for ee=1:mesh{lev}.NE
%            val=is_inside_the_square(bounding_box{pp},mean(mesh{lev}.node(mesh{lev}.edge(ee,:),:),1)) ;
%          end
         
         for nn=1:mesh{lev}.N
           if(is_inside_the_square(bounding_box{pp},mesh{lev}.node(nn,:) ) )
               vertices_per_processor{lev}{pp}=[vertices_per_processor{lev}{pp};nn];
           end
         end         

     end
 end
 
maps.Patch_Node{lev}{nn}
 
end

