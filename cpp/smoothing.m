close all
smoothing_steps=1000;
num_procs=4;
num_colors=7;
clear V;
clear male;
clear residual;
clear energy;
clear error;
A=Amat;
b=resfile;
sol=A\b;
for cc=1:num_colors
    for pp=1:num_procs   
           vertices_cc=color2vertex{pp}{cc};
           pp_dofs{pp}=[];
           for ii=1:length(vertices_cc)
               pp_dofs{pp}=[pp_dofs{pp},Patch{pp}{vertices_cc(ii)}];               
           end
           pp_dofs{pp}=unique(pp_dofs{pp});
    end
     for pp1=1:num_procs   
           for pp2=1:num_procs   
               V{cc}(pp1,pp2)=length(intersect(pp_dofs{pp1},pp_dofs{pp2}));
            end
    end   
    
end



for cc=1:num_colors
cont=0;
     Matr=zeros();
     for pp1=1:num_procs   
           vertices_pp1=color2vertex{pp1}{cc};
           for ii1=1:length(vertices_pp1)
           patch1=Patch{pp1}{vertices_pp1(ii1)};
           for pp2=(pp1+1):num_procs   
                 inter=intersect(pp_dofs{pp1},pp_dofs{pp2});
                 vertices_pp2=color2vertex{pp2}{cc};
                 for ii2=1:length(vertices_pp2)
                     patch2=Patch{pp2}{vertices_pp2(ii2)};
                     if(length(intersect(patch1,patch2))>0)
                         cont=cont+1;
                         
                         male{cc,cont}=[pp1,pp2,vertices_pp1(ii1),vertices_pp2(ii2)];
                     end
                 end
           end
           
          end
    end   
    
end






kk=0;
x=zeros(length(A(:,1)),1);
res=b;
% smoothing steps
for ss=1:smoothing_steps
    % loop on all the colors
    for cc=1:num_colors
        
    % loop on all the processors (simulate parallelization)
       for pp=1:num_procs
           correction{pp}=zeros(length(A(:,1)),1);    
           vertices_cc=color2vertex{pp}{cc};
           for ii=1:length(vertices_cc)
               patch_dofs=Patch{pp}{vertices_cc(ii)};
               b_Loc=res(patch_dofs)-A(patch_dofs,:)*correction{pp};
               x_loc=A(patch_dofs,patch_dofs)\b_Loc;
               correction{pp}(patch_dofs)=correction{pp}(patch_dofs)+x_loc;
           end
           %x=x+correction{pp};
       end
       
       for pp=1:num_procs
          x=x+correction{pp};
       end
       
       sum_norm_corr=norm(correction{1},1)+norm(correction{2},1)+norm(correction{3},1)+norm(correction{4},1);
       norm_summ_corr=norm((correction{1}+correction{2}+correction{3}+correction{4}),1);
       if(abs(norm_summ_corr-sum_norm_corr)>0.000001)
           fermami=1
       end
  
       kk=kk+1;
       res=b-A*x;
       residual(kk)=norm(res);
       energy(kk)=0.5*x'*A*x-b'*x;
       if(kk>1)
           if(energy(kk)<=energy(kk-1))
             energy_decreasing(kk-1)=1;
           else
             energy_decreasing(kk-1)=0;
           end
       end
       error(kk)=norm(x-sol);
       DotsProduct=zeros(num_procs,num_procs);
       for pp1=1:num_procs
           for pp2=1:num_procs
               DotsProduct(pp1,pp2)=correction{pp1}'*correction{pp2};
           end
       end
       
    end
    
    
    
end

plot(log10(residual))
figure
plot((energy))
figure
plot(log10(error))
figure
plot((energy_decreasing))