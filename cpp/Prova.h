#include <dolfin.h>

using namespace dolfin;



 void coloring(std::map<int, std::set<unsigned int> >&  shared_vertices,  std::shared_ptr<dolfin::Mesh> mesh, std::vector<std::vector< int > >& topology_N2N,
              std::vector<unsigned int>& vertex_color,std::vector<unsigned int>& vertex_global_dof,
              std::vector<unsigned int>& vertex_shared_color,std::vector<unsigned int>& vertex_shared_global_dof,
              std::map<unsigned int, unsigned int>& shared_rank_map, std::vector<std::vector<unsigned int> >& vertex_shared_global_dof_and_color,
              std::vector<unsigned int>& used_color)
 {



unsigned int count_vertex=0;
int world_rank,world_size;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
// std::vector<std::set<unsigned int> > vals;
 //              shared_rank.push_back(*it_procs); //contains all the communicating processes (multiple times)
      
// std::vector<unsigned int> vertex_shared_color;
// std::vector<unsigned int> vertex_shared_global_dof;
// std::vector<unsigned int> vertex_color(mesh->num_vertices(),0);
// std::vector<unsigned int> vertex_global_dof(mesh->num_vertices(),0);

std::cout<<"world_rank=="<<world_rank<<", PRE coloring"<<std::endl;

   for (std::map< int, std::set<unsigned int> >::iterator shared_node=shared_vertices.begin(); shared_node!= shared_vertices.end(); ++shared_node)
     { 
     
        auto shared_procs=shared_node->second;
        std::vector<unsigned int> patch_shared_color;
        unsigned int vertex_index=shared_node->first;
        auto actual_vertex=Vertex(*mesh, vertex_index);
        unsigned int global_vertex_index=actual_vertex.global_index();
          // for the first node, use the color=1
          //std::cout<<"world_rank=="<<world_rank<<", vertex_index=="<<vertex_index<<", global_vertex_index "<< global_vertex_index<<std::endl;
        if(vertex_color[vertex_index]==0)
         {
         
         
         if(used_color.size()==1)
           {
            used_color.push_back(1);
            vertex_color[vertex_index]=1;
            vertex_global_dof[vertex_index]=global_vertex_index;
            vertex_shared_color.push_back(1);
            vertex_shared_global_dof.push_back(global_vertex_index);
            for ( std::set<unsigned int>::iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
                 {vertex_shared_global_dof_and_color[shared_rank_map[*it_procs]].push_back(global_vertex_index);
                  vertex_shared_global_dof_and_color[shared_rank_map[*it_procs]].push_back(vertex_color[vertex_index]);
                 }           
            
            count_vertex++;
           }
         else
         {
      	   auto N2N=topology_N2N[vertex_index];
     	   
    	    // check which color I can use 
    	    for(int ii=0;ii<N2N.size();ii++)
   	        {
   	          // consider all the colors >1 of the patch
  	          if(vertex_color[N2N[ii]]>0)
  	            patch_shared_color.push_back(vertex_color[N2N[ii]]);
  	         }
  	         
  	     // sort and unique patch_shared_color
 		 std::sort(patch_shared_color.begin(),patch_shared_color.end());
         auto patch_shared_color_tmp = std::unique(patch_shared_color.begin(), patch_shared_color.end());
    	 patch_shared_color.erase(patch_shared_color_tmp, patch_shared_color.end());   	
     	 
    	 // if the numer of colors up to now used (except the zero) is equal to the ones on the patch
    	 // then add onother color
    	 unsigned int used_color_size=used_color.size();
    	 if(patch_shared_color.size()==used_color_size-1)   
    	   { 
    	    used_color.push_back(1);
    	    vertex_color[vertex_index]=used_color_size;
    	    
    	    vertex_global_dof[vertex_index]=global_vertex_index;
    	    vertex_shared_color.push_back(used_color_size);
            vertex_shared_global_dof.push_back(global_vertex_index);
            for ( std::set<unsigned int>::iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
                 {vertex_shared_global_dof_and_color[shared_rank_map[*it_procs]].push_back(global_vertex_index);
                  vertex_shared_global_dof_and_color[shared_rank_map[*it_procs]].push_back(vertex_color[vertex_index]);
                 }
            count_vertex++; 
    	   }     
    	 //otherwise we can opt among one of the colors already use
    	 // we discard the ones in patch_shared_color
    	 // and of the remaining in used_colors, we take the one which has less vertices
    	 else
    	   {
    	    std::vector<unsigned int> unused_colors(used_color);
    	    std::vector<int> unused_colors_range(unused_colors.size());
            std::iota(unused_colors_range.begin(), unused_colors_range.end(), 0);
            
    	    for(int ii=0;ii<patch_shared_color.size();ii++)
                {unused_colors.erase( unused_colors.begin() +  patch_shared_color[ii]-ii);
                 unused_colors_range.erase( unused_colors_range.begin() +  patch_shared_color[ii]-ii);
                }
            vertex_color[vertex_index]=unused_colors_range[std::distance(unused_colors.begin(), std::min_element(unused_colors.begin(), unused_colors.end()))];
            vertex_global_dof[vertex_index]=global_vertex_index;
    	    vertex_shared_color.push_back(vertex_color[vertex_index]);
            vertex_shared_global_dof.push_back(global_vertex_index);
            count_vertex++; 
            for ( std::set<unsigned int>::iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
                 {vertex_shared_global_dof_and_color[shared_rank_map[*it_procs]].push_back(global_vertex_index);
                  vertex_shared_global_dof_and_color[shared_rank_map[*it_procs]].push_back(vertex_color[vertex_index]);
                 }
            }
                                 
          }
          }
         patch_shared_color.clear();
     }    
std::cout<<"world_rank=="<<world_rank<<", POST coloring"<<std::endl;

}   




std::vector<std::vector< int > > topologyN2N( std::shared_ptr<dolfin::Mesh> mesh, dolfin::MeshConnectivity &topology_N2F, dolfin::MeshConnectivity &topology_F2N)
{

 	std::vector<std::vector< int > > topology_N2N(mesh->num_vertices());
	unsigned int gdim = mesh->geometry().dim();  
    int world_rank,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	for(int rr=0;rr<world_size;rr++)
	{
	if(world_rank==rr)
	{
	for(int nn=0;nn<mesh->num_vertices();nn++)
	{
	auto actual_vertex=Vertex(*mesh, nn);
	auto point_vertex=actual_vertex.point();
	//std::cout<<"world_rank=="<<world_rank<<", point_vertex==("<<point_vertex[0]<<", "<<point_vertex[1]<<")"<<std::endl;


	  for (MeshEntityIterator ee(  MeshEntity(*mesh, 0, nn),gdim-1 ); !ee.end(); ++ee)
	 {
		auto edge_dof=topology_N2F(nn)[ee.pos()];
		for(int nn2=0;nn2<gdim;nn2++)
		   if(nn!=topology_F2N(edge_dof)[nn2])
		   {
			topology_N2N[nn].push_back(topology_F2N(edge_dof)[nn2]);
			auto patch_vertex=Vertex(*mesh, topology_F2N(edge_dof)[nn2]);
			auto point_patch=patch_vertex.point();
		   // std::cout<<"world_rank=="<<world_rank<<", point_patch==("<<point_patch[0]<<", "<<point_patch[1]<<")"<<std::endl;

		  }            
	 }
	}

	}
	MPI_Barrier(MPI_COMM_WORLD);
	}  
	
	return topology_N2N;
}   


std::vector<std::vector< int > > topologyN2PatchDofs( std::shared_ptr<dolfin::Mesh> mesh,std::shared_ptr<const dolfin::GenericDofMap> dofmap, dolfin::MeshConnectivity &topology_N2F)
{
	std::vector<std::vector< int > > Patch(mesh->num_vertices());
    unsigned int gdim = mesh->geometry().dim();  
    VertexIterator vertex=VertexIterator(*mesh);
    
	for (; !vertex.end(); ++vertex)
	{
	auto vertex_index=vertex->index();
	std::vector<long unsigned int> vertex_vector(1);
	vertex_vector[0]=vertex_index;
	auto actual_vertex=Vertex(*mesh, vertex_index);
	auto point_vertex=actual_vertex.point();

	// add to the patch the dofs related to the node
	auto tmp_node=dofmap->entity_dofs(*mesh, 0,vertex_vector);
	for(int ii=0;ii<tmp_node.size();ii++)
		Patch[vertex_index].push_back(tmp_node[ii]);

	for (MeshEntityIterator ee(  MeshEntity(*mesh, 0, vertex_index),gdim-1 ); !ee.end(); ++ee)
	{
	  auto edge_dof=topology_N2F(vertex_index)[ee.pos()];
	  std::vector<long unsigned int> edge_vector(1);
	  edge_vector[0]=edge_dof;
	  auto actual_edge=Edge(*mesh,edge_dof);
	  auto point_edge=actual_edge.midpoint();
	  // add to the patch the dofs related to the faces connected to the node
	  auto tmp_face=dofmap->entity_dofs(*mesh, 1,edge_vector);
	  for(int ii=0;ii<tmp_face.size();ii++)
		  Patch[vertex_index].push_back(tmp_face[ii]);
	  //std::cout<<"NOT SHARED world_rank: "<<world_rank<<", topology_N2F: " <<vertex_index<<" coord "<<point_vertex[0]<<", "<<point_vertex[1]<< ", edge: "<<edge_dof<<" coord "<< point_edge[0]<<", "<<point_edge[1] <<std::endl;
	}
	}
	
	
	// loop on all the vertices (also ghost)
// while(it_sharedvertex != shared_vertices.end())
// {
// auto vertex_index=it_sharedvertex->first;
// auto actual_vertex=Vertex(*mesh, vertex_index);
// auto point_vertex=actual_vertex.point();
// 
// for (MeshEntityIterator ee(  MeshEntity(*mesh, 0, it_sharedvertex->first),gdim-1 ); !ee.end(); ++ee)
// {
//   auto edge_num=topology_N2F(it_sharedvertex->first)[ee.pos()];
//   auto edge=Edge(*mesh,edge_num);
//   
//   auto edge_dof=topology_N2F(vertex_index)[ee.pos()];
//   auto actual_edge=Edge(*mesh,edge_dof);
//   auto point_edge=actual_edge.midpoint();
//  // std::cout<<"SHARED world_rank: "<<world_rank<<", topology_N2F: " <<it_sharedvertex->first<<" coord "<<point_vertex[0]<<", "<<point_vertex[1]<< ", edge: "<<edge_dof<<" coord "<< point_edge[0]<<", "<<point_edge[1] <<std::endl;
// }
// it_sharedvertex++;
// }

return Patch;
}










// loop on all the vertices (not ghost, but the patch will consider also the ghost faces)
// for (; !vertex.end(); ++vertex)
// {
// auto vertex_index=vertex->index();
// std::vector<long unsigned int> vertex_vector(1);
// vertex_vector[0]=vertex_index;
// auto actual_vertex=Vertex(*mesh, vertex_index);
// auto point_vertex=actual_vertex.point();
// 
// // add to the patch the dofs related to the node
// auto tmp_node=dofmap->entity_dofs(*mesh, 0,vertex_vector);
// for(int ii=0;ii<tmp_node.size();ii++)
//     Patch[vertex_index].push_back(tmp_node[ii]);
// 
// for (MeshEntityIterator ee(  MeshEntity(*mesh, 0, vertex_index),gdim-1 ); !ee.end(); ++ee)
// {
//   auto edge_dof=topology_N2F(vertex_index)[ee.pos()];
//   std::vector<long unsigned int> edge_vector(1);
//   edge_vector[0]=edge_dof;
//   auto actual_edge=Edge(*mesh,edge_dof);
//   auto point_edge=actual_edge.midpoint();
//   // add to the patch the dofs related to the faces connected to the node
//   auto tmp_face=dofmap->entity_dofs(*mesh, 1,edge_vector);
//   for(int ii=0;ii<tmp_face.size();ii++)
//       Patch[vertex_index].push_back(tmp_face[ii]);
//   //std::cout<<"NOT SHARED world_rank: "<<world_rank<<", topology_N2F: " <<vertex_index<<" coord "<<point_vertex[0]<<", "<<point_vertex[1]<< ", edge: "<<edge_dof<<" coord "<< point_edge[0]<<", "<<point_edge[1] <<std::endl;
// }
// }

// loop on all the vertices (also ghost)
// while(it_sharedvertex != shared_vertices.end())
// {
// auto vertex_index=it_sharedvertex->first;
// auto actual_vertex=Vertex(*mesh, vertex_index);
// auto point_vertex=actual_vertex.point();
// 
// for (MeshEntityIterator ee(  MeshEntity(*mesh, 0, it_sharedvertex->first),gdim-1 ); !ee.end(); ++ee)
// {
//   auto edge_num=topology_N2F(it_sharedvertex->first)[ee.pos()];
//   auto edge=Edge(*mesh,edge_num);
//   
//   auto edge_dof=topology_N2F(vertex_index)[ee.pos()];
//   auto actual_edge=Edge(*mesh,edge_dof);
//   auto point_edge=actual_edge.midpoint();
//  // std::cout<<"SHARED world_rank: "<<world_rank<<", topology_N2F: " <<it_sharedvertex->first<<" coord "<<point_vertex[0]<<", "<<point_vertex[1]<< ", edge: "<<edge_dof<<" coord "<< point_edge[0]<<", "<<point_edge[1] <<std::endl;
// }
// it_sharedvertex++;
// }


// loop on the vertices
// for(int nn=0;nn<mesh->num_vertices() ;nn++)
//     {// loop on the dofs of the given vertex
//     auto actual_vertex=Vertex(*mesh, nn);
//     auto point_vertex=actual_vertex.point();
//     for(int ii=0;ii<Patch[nn].size();ii++)
//         {
//         // std::cout<<"NOT SHARED["<<world_rank<<"], vertex="<<nn<<",with coord["<< point_vertex[0]<<", "<<point_vertex[1]<<"], and dof["<<ii<<"]="<<Patch[nn][ii]<<std::endl;
//         }
//     }