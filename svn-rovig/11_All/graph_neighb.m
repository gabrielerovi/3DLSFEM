function graph=graph_neighb(mesh)

N=mesh.N;
T=mesh.NT;
node_per_elem=mesh.node_per_elem;
node=mesh.node;
N_to_T=mesh.N_to_T;
cont=0;
vertex=1;
flag_node(1)=1;
flag_node=zeros(N,1);
flag_elem=zeros(T,1);


while(cont<N)

    % we define a temporary element flag vector
    tmp_flag_elem=zeros(T,1);
    
    % here we fix a patch and consider all its non-marked nodes that will
    % be added to the graph
    % at the next while-iteration, we need a new node. We are going to
    % consider one belonging to the boundary of the patch such that it has
    % at least one neighbour element with at least one non-marked node.
    % If this does not happen, search for a node that is non marked and
    % restart
    
    elements=N_to_T{vertex};
    
    for ee=1 : length(elements)
        
        elem=elements{ee};
        
        %let us temporarly mark this element for later use
        tmp_flag_elem(elem)=1;
        % if the element has not been marked yet
        if( flag_elem(elem)==0)
            % then let us consider this element elem and mark it
            flag_elem(elem)=1;
           
            elem_nodes=mesh.elem(elem,:);
          % consider the nodes of the element
          for ii=1:node_per_elem
              
              % if the node has not been marked yet
              if(flag_node(elem_nodes(ii))==0)
                 % mark it
                 flag_node(elem_nodes(ii))=1;
                 % and add to the graph the new node
                 cont=cont+1;
                 graph(cont)=elem_nodes(ii);
              end
              
              
          end      
            
        end
        
    end
    

        please_exit_me=false;
        elements=N_to_T{vertex}; 
    
        % search again among the elements of the patch
        for ee=1 : length(elements)
        elem_nodes=mesh.elem(elem,:);
        % for a fixed element, consider all its nodes
          for ii=1:node_per_elem  
              
              %then consider its neighbour elements but discard the
              %elements previously examined
              tmp_elements=N_to_T{elem_nodes(ii)};
              
              for tmp_ee=1 : length(tmp_elements)
                  
                  % discard the elements previously examined
                  if(tmp_flag_elem(tmp_elements{tmp_ee})==0)
                      
                      % if the element has not been examined yet, maybe
                      % also its nodes are not
                      tmp_nodes=mesh.elem(tmp_elements{tmp_ee},:); 
                      
                      for tmp_nn=1:node_per_elem
                          if(flag_node(tmp_nodes(tmp_nn))==0)
                          cont=cont+1;
                          vertex=tmp_nodes(tmp_nn);
                          flag_node(vertex)=1;
                          graph(cont)=vertex;
                          please_exit_me=true;
                          end
                          if(please_exit_me==true)
                              break;
                          end
                           if(please_exit_me==true)
                              break;
                           end
                      end
                      
                  end
                  
                  
              if(please_exit_me==true)
                break;
              end  
              
              end
              
           if(please_exit_me==true)
            break;
           end
          
          end
          
          
          if(please_exit_me==true)
            break;
          end
          
        end
        
        if(please_exit_me==false && cont<N)
            node_found=false;
            cont1=0;
           while(node_found==false)
               cont1=cont1+1;
               if(flag_node(cont1)==0)
                   node_found=true;
                   vertex=cont1;
               end
           end
        end
        
end
              
              
              
              
              
              
              
              
              
              
              
      
    
end