function [edgegraphinterior,edgegraphgammac]=graph_neighb_interior_and_gammac(mesh)

N=mesh.N;
NC=length(mesh.N_contact);
T=mesh.NT;
node_per_elem=mesh.node_per_elem;
node=mesh.node;
N_to_T=mesh.N_to_T;
cont=0;
vertex=1;
flag_node(1)=1;
flag_node=zeros(N,1);
flag_elem=zeros(T,1);


% we do not consider nodes on gammaC
flag_node(mesh.N_contact)=1;

cont=0;

N_dirichlet=mesh.N_dirichlet;

NE=length(mesh.edge);
N_bool_contact=zeros(mesh.N,1);
N_bool_contact(mesh.N_contact)=1;
E_bool_contact=zeros(NE,1);

NEContact=0;
flag_edge=zeros(NE,1);

E_to_T=cell(NE,1);
cont=0;
for ee=1:NE
edge=mesh.edge(ee,:);

if(N_bool_contact(edge(1))==1|| N_bool_contact(edge(2))==1)
E_bool_contact=1;
NEContact=NEContact+1;
flag_edge(ee)=1;
cont=cont+1;
E_contact(cont)=ee;
end

E_to_T{ee}=intersect(cell2mat(mesh.N_to_T{edge(1)}),cell2mat(mesh.N_to_T{edge(2)}));

end






list_of_edges_added_now=[];
edge=1;
cont=0;
 while(cont<NE-NEContact)

       
       T=E_to_T{edge};
       list_of_nodes_added_now=[];
       for tt=1:length(T)
       elem=T(tt); 
       elem_edges=mesh.elemE(elem,:);
       
       for ii=1:length(elem_edges)
           if(flag_edge(elem_edges(ii))==0)
             cont=cont+1;
             flag_edge(elem_edges(ii))=1;
             edgegraphinterior(cont)=elem_edges(ii);
             list_of_edges_added_now=[list_of_edges_added_now;elem_edges(ii)];
           end
       end
       end
       
       
        please_exit=false;      
        for jj=1:length(list_of_edges_added_now)     
          T=E_to_T(list_of_edges_added_now(jj));  
          for tt=1:length(T)
           elem=T{tt}; 
           elem_edges=mesh.elemE(elem,:);              
              for ii=1:length(elem_edges)
                       if(flag_edge(elem_edges(ii))==0)
                           edge=elem_edges(ii);
                           please_exit=true; 
                       end
                       if(please_exit==true)
                           break;
                       end
              end
              
              
              if(please_exit==true)
                break;
              end
                       
          end
          
          
          
          if(please_exit==true)
           break;
          end
                       
                       
        end
        
        if(please_exit==false)
          edge=min(find(flag_edge==0))  ;
        end
        
  end








flag_edge=ones(NE,1);
flag_edge(E_contact)=0;


% we only consider nodes on gammaC

edge=(1);
cont=0;
  while(cont<NEContact)

       
       T=E_to_T{edge};
       list_of_edges_added_now=[];
       for tt=1:length(T)
       elem=T(tt); 
       elem_edges=mesh.elemE(elem,:);
       
       for ii=1:length(elem_edges)
           if(flag_edge(elem_edges(ii))==0)
             cont=cont+1;
             flag_edge(elem_edges(ii))=1;
             edgegraphgammac(cont)=elem_edges(ii);
             list_of_edges_added_now=[list_of_edges_added_now;elem_edges(ii)];
           end
       end
       end
       
       
        please_exit=false;      
        for jj=1:length(list_of_edges_added_now)     
          T=E_to_T{list_of_edges_added_now(jj)};  
          for tt=1:length(T)
           elem=T(tt); 
           elem_edges=mesh.elemE(elem,:);              
              for ii=1:length(elem_edges)
                       if(flag_edge(elem_edges(ii))==0)
                           edge=elem_edges(ii);
                           please_exit=true; 
                       end
                       if(please_exit==true)
                           break;
                       end
              end
              
              
              if(please_exit==true)
                break;
              end
                       
          end
          
          
          
          if(please_exit==true)
           break;
          end
                       
                       
        end
        
        if(please_exit==false)
          edge=min(find(flag_edge==0))  ;
        end
        
  end
  
  
  
  
  
  
  
  































% while(cont<N-NC)
% 
%     % we define a temporary element flag vector
%     tmp_flag_elem=zeros(T,1);
%     
%     % here we fix a patch and consider all its non-marked nodes that will
%     % be added to the graphinterior
%     % at the next while-iteration, we need a new node. We are going to
%     % consider one belonging to the boundary of the patch such that it has
%     % at least one neighbour element with at least one non-marked node.
%     % If this does not happen, search for a node that is non marked and
%     % restart
%     
%     elements=N_to_T{vertex};
%     
%     for ee=1 : length(elements)
%         
%         elem=elements{ee};
%         
%         %let us temporarly mark this element for later use
%         tmp_flag_elem(elem)=1;
%         % if the element has not been marked yet
%         if( flag_elem(elem)==0)
%             % then let us consider this element elem and mark it
%             flag_elem(elem)=1;
%            
%             elem_nodes=mesh.elem(elem,:);
%           % consider the nodes of the element
%           for ii=1:node_per_elem
%               
%               % if the node has not been marked yet
%               if(flag_node(elem_nodes(ii))==0)
%                  % mark it
%                  flag_node(elem_nodes(ii))=1;
%                  % and add to the graphinterior the new node
%                  cont=cont+1;
%                  graph.graphinterior(cont)=elem_nodes(ii);
%               end
%               
%               
%           end      
%             
%         end
%         
%     end
%     
% 
%         please_exit_me=false;
%         elements=N_to_T{vertex}; 
%     
%         % search again among the elements of the patch
%         for ee=1 : length(elements)
%         elem_nodes=mesh.elem(elem,:);
%         % for a fixed element, consider all its nodes
%           for ii=1:node_per_elem  
%               
%               %then consider its neighbour elements but discard the
%               %elements previously examined
%               tmp_elements=N_to_T{elem_nodes(ii)};
%               
%               for tmp_ee=1 : length(tmp_elements)
%                   
%                   % discard the elements previously examined
%                   if(tmp_flag_elem(tmp_elements{tmp_ee})==0)
%                       
%                       % if the element has not been examined yet, maybe
%                       % also its nodes are not
%                       tmp_nodes=mesh.elem(tmp_elements{tmp_ee},:); 
%                       
%                       for tmp_nn=1:node_per_elem
%                           if(flag_node(tmp_nodes(tmp_nn))==0)
%                           cont=cont+1;
%                           vertex=tmp_nodes(tmp_nn);
%                           flag_node(vertex)=1;
%                           graph.graphinterior(cont)=vertex;
%                           please_exit_me=true;
%                           end
%                           if(please_exit_me==true)
%                               break;
%                           end
%                            if(please_exit_me==true)
%                               break;
%                            end
%                       end
%                       
%                   end
%                   
%                   
%               if(please_exit_me==true)
%                 break;
%               end  
%               
%               end
%               
%            if(please_exit_me==true)
%             break;
%            end
%           
%           end
%           
%           
%           if(please_exit_me==true)
%             break;
%           end
%           
%         end
%         
%         if(please_exit_me==false && cont<N-NC)
%             node_found=false;
%             cont1=0;
%            while(node_found==false)
%                cont1=cont1+1
%                if(flag_node(cont1)==0)
%                    node_found=true;
%                    vertex=cont1;
%                end
%            end
%         end
%         
% end
%               
%               
%               
%               
% flag_node=ones(N,1);
% flag_elem=zeros(T,1);
% 
% 
% % we only consider nodes on gammaC
% N_contact=mesh.N_contact;
% flag_node(N_contact)=0;  
% 
% vertex=N_contact(1);
% cont=0;
%   while(cont<NC)
% 
%        N_contact(1)
%        
%        T=N_to_T{vertex};
%        list_of_nodes_added_now=[];
%        for tt=1:length(T)
%        elem=T{tt}; 
%        elem_nodes=mesh.elem(elem,:);
%        
%        for ii=1:length(elem_nodes)
%            if(flag_node(elem_nodes(ii))==0)
%              cont=cont+1;
%              flag_node(elem_nodes(ii))=1;
%              graph.graphgammac(cont)=elem_nodes(ii);
%              list_of_nodes_added_now=[list_of_nodes_added_now;elem_nodes(ii)];
%            end
%        end
%        end
%        
%        
%         please_exit=false;      
%         for jj=1:length(list_of_nodes_added_now)     
%           T=N_to_T{list_of_nodes_added_now(jj)};  
%           for tt=1:length(T)
%            elem=T{tt}; 
%            elem_nodes=mesh.elem(elem,:);              
%               for ii=1:length(elem_nodes)
%                        if(flag_node(elem_nodes(ii))==0)
%                            vertex=elem_nodes(ii);
%                            please_exit=true; 
%                        end
%                        if(please_exit==true)
%                            break;
%                        end
%               end
%               
%               
%               if(please_exit==true)
%                 break;
%               end
%                        
%           end
%           
%           
%           
%           if(please_exit==true)
%            break;
%           end
%                        
%                        
%         end
%         
%         if(please_exit==false)
%           vertex=min(find(flag_node==0))  ;
%         end
%         
%   end
        
              
              
              
      
    
end