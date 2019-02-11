#include <dolfin.h>
#include <fstream> 
#include "writetofile.hpp"
using namespace dolfin;

 void coloring(const std::map<int, std::set<unsigned int> >& shared_vertices,  std::shared_ptr<dolfin::Mesh> mesh, 
              const MeshConnectivity & topology_N2N,std::vector<unsigned int>& vertex_color,
              //std::vector<unsigned int>& vertex_global_dof,std::vector<unsigned int>& vertex_shared_color,
              std::vector<unsigned int>& vertex_shared_global_dof,const std::map<unsigned int, unsigned int>& shared_rank_map, 
              std::vector<std::vector<unsigned int> >& vertex_shared_global_dof_and_color,std::vector<unsigned int>& used_color)
             // std::vector<bool> &use_vertex_color)
 {



unsigned int count_vertex=0;
int world_rank,world_size;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
MPI_Comm_size(MPI_COMM_WORLD, &world_size);

   // loop on all the shared nodes
   for (std::map< int, std::set<unsigned int> >::const_iterator shared_node=shared_vertices.begin(); shared_node!= shared_vertices.end(); ++shared_node)
     { 
     
        auto shared_procs=shared_node->second;
        std::vector<unsigned int> patch_shared_color;
        unsigned int vertex_index=shared_node->first;
        auto actual_vertex=Vertex(*mesh, vertex_index);
        unsigned int global_vertex_index=actual_vertex.global_index();
        // no previous process has colored the actual process, therefore we can add it 
        if(vertex_color[vertex_index]==0)
         {
         
         
         if(used_color.size()==1)
           {
            used_color.push_back(1);
            vertex_color[vertex_index]=1;
           // vertex_global_dof[vertex_index]=global_vertex_index;
            //vertex_shared_color.push_back(1);
            vertex_shared_global_dof.push_back(global_vertex_index);
            for ( std::set<unsigned int>::const_iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
                 {vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(global_vertex_index);
                  vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(vertex_color[vertex_index]);
                 }           
            
            count_vertex++;
           }
         else
         {           
      	   auto N2N=topology_N2N(vertex_index);
     	   // check which color I can use 
     	   for (MeshEntityIterator edge(  MeshEntity(*mesh, 0, vertex_index),1 ); !edge.end(); ++edge)
	            if(vertex_color[N2N[edge.pos()]]>0)
  	               patch_shared_color.push_back(vertex_color[N2N[edge.pos()]]);

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
    	    
    	   // vertex_global_dof[vertex_index]=global_vertex_index;
    	   // vertex_shared_color.push_back(used_color_size);
            vertex_shared_global_dof.push_back(global_vertex_index);
            for ( std::set<unsigned int>::const_iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
                 {vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(global_vertex_index);
                  vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(vertex_color[vertex_index]);
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
          //  vertex_global_dof[vertex_index]=global_vertex_index;
    	 //   vertex_shared_color.push_back(vertex_color[vertex_index]);
            vertex_shared_global_dof.push_back(global_vertex_index);
            count_vertex++; 
            for ( std::set<unsigned int>::const_iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
                 {vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(global_vertex_index);
                  vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(vertex_color[vertex_index]);
                 }
            }
                                 
          }
          }

         patch_shared_color.clear();
     }    

}   
























std::map<unsigned int, unsigned int> global2local_vertex(std::shared_ptr<dolfin::Mesh> mesh, const std::map<int, std::set<unsigned int> >& shared_vertices)
{ 
std::map<unsigned int, unsigned int> global_to_local_vertex;   
 for (std::map< int, std::set<unsigned int> >::const_iterator shared_node=shared_vertices.begin(); shared_node!= shared_vertices.end(); ++shared_node)
 {
   auto vertex_index=shared_node->first;
   auto actual_vertex=Vertex(*mesh, vertex_index);
   global_to_local_vertex[actual_vertex.global_index()] = vertex_index;
    
 }
 return global_to_local_vertex;
}













void find_communicating_processes(std::vector<unsigned int> &shared_rank, 
                                  std::vector<unsigned int> &all_shared_ranks, 
                                  const std::map<int, std::set<unsigned int> >& shared_vertices, 
                                  std::map<unsigned int, unsigned int> &shared_rank_map)
{
int world_rank,world_size;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
MPI_Comm_size(MPI_COMM_WORLD, &world_size);

//std::vector<unsigned int> keys; 
std::vector<std::set<unsigned int> > vals;

unsigned int cont=0;
     
    for (std::map< int, std::set<unsigned int> >::const_iterator it=shared_vertices.begin(); it!= shared_vertices.end(); ++it)
      {
       //keys.push_back(it->first);
       vals.push_back(it->second);
       for ( std::set<unsigned int>::iterator setit = vals[cont].begin(); setit != vals[cont].end(); ++setit)
             shared_rank.push_back(*setit); //contains all the communicating processes (multiple times)
             
       cont++;
       }

// now we sort and make unique shared_rank, which contains only the communicating processes
// while all_shared_ranks contains also the current processor
std::sort(shared_rank.begin(),shared_rank.end());
auto shared_rank_tmp = std::unique(shared_rank.begin(), shared_rank.end());
shared_rank.erase(shared_rank_tmp, shared_rank.end()); 
all_shared_ranks=shared_rank;
all_shared_ranks.push_back(world_rank);
std::sort(all_shared_ranks.begin(),all_shared_ranks.end());


for(int ii=0;ii<shared_rank.size();ii++)
   shared_rank_map[shared_rank[ii]]=ii;
}












std::vector<bool>  divide_nodes_among_processes(std::shared_ptr<dolfin::Mesh> mesh,const std::map<int, std::set<unsigned int> >& shared_vertices,
                                  const std::vector<unsigned int> &shared_rank,const std::vector<unsigned int> &all_shared_ranks, 
                                  const std::map<unsigned int, unsigned int> &shared_rank_map, const std::map<unsigned int, unsigned int> &global_to_local_vertex)
{
int world_rank,world_size;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 
MPI_Comm_size(MPI_COMM_WORLD, &world_size);  
unsigned int num_domain_vertices=mesh->topology().ghost_offset(0);
unsigned int num_all_vertices=mesh->num_vertices();
std::vector<bool> use_vertex_color(num_all_vertices,false);

  std::vector<std::vector<unsigned int> > buffer_vector(shared_rank.size());
    
    for (std::map< int, std::set<unsigned int> >::const_iterator vertex_it=shared_vertices.begin(); vertex_it!= shared_vertices.end(); ++vertex_it)
      {
       unsigned int vertex_index=vertex_it->first;
       auto actual_vertex=Vertex(*mesh, vertex_index);
       unsigned int global_vertex_index=actual_vertex.global_index();
       auto vertex_processes=vertex_it->second;
       for ( std::set<unsigned int>::iterator comm_procs_it = vertex_processes.begin(); comm_procs_it != vertex_processes.end(); ++comm_procs_it)
             {
             unsigned int rank_to_send=*comm_procs_it;
             buffer_vector[shared_rank_map.at(rank_to_send)].push_back(global_vertex_index);
             // 0 if it belongs to world_rank domain AND if world_rank procs is smaller than the one receiving
             // in this case, the receving proc will not consider this node as its
             if(vertex_index<num_domain_vertices && world_rank<rank_to_send)
             buffer_vector[shared_rank_map.at(rank_to_send)].push_back( 0 ); 
             // otherwise if it outside the domain or if the receiving procs has a lower rank than world_rank, it will own this node 
             else
             buffer_vector[shared_rank_map.at(rank_to_send)].push_back( 1 ); // it does not belong to the other process
             }

       }

for(int ii=0;ii<num_domain_vertices;ii++)
use_vertex_color[ii]=true;
    
    for(int ii=0;ii<shared_rank.size();ii++)
    {
        {
        MPI_Request request;
        unsigned int* buffer;
            int tag=world_rank*world_size+shared_rank[ii];
            int rank2send=shared_rank_map.at(shared_rank[ii]);
            int send_rank=shared_rank[ii];
            int send_size=buffer_vector[ii].size();
            ///std::cout<<"send rank=="<<world_rank<<", to="<<shared_rank[ii]<<std::endl;
            buffer=(unsigned int*)malloc(sizeof(unsigned int)*buffer_vector[ii].size());
            for(int jj=0;jj<send_size;jj++)
                buffer[jj]=buffer_vector[ii][jj];
            MPI_Isend(buffer,send_size,MPI_UNSIGNED,send_rank,tag,MPI_COMM_WORLD,&request);
        }
    }

 // you must do a scatter of all num_domain_vertices
    for(int ii=0;ii<shared_rank.size();ii++)
    {// receive from processes with a lower rank
        MPI_Status status_now;
            int recv_size;
            unsigned int* buffer_recv;
            //std::cout<<"recv rank=="<<world_rank<<", from="<<shared_rank[ii]<<std::endl;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status_now);
            MPI_Get_count(&status_now, MPI_UNSIGNED, &recv_size);
            buffer_recv=(unsigned int*)malloc(sizeof(unsigned int)*recv_size);
            MPI_Recv(buffer_recv, recv_size, MPI_UNSIGNED, status_now.MPI_SOURCE, status_now.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);// &status_now);
            
            for(int jj=0;jj<recv_size;jj=jj+2)    
            {
            
            unsigned int local_index=global_to_local_vertex.at(buffer_recv[jj]);
            if( (use_vertex_color[local_index]==true&&buffer_recv[jj+1]==0) || local_index>=num_domain_vertices)
            use_vertex_color[local_index]=false;
            //if(world_rank==1)
            //   std::cout<<"-------------world_rank="<<world_rank<<", recv_size="<<recv_size<<", use_vertex_color["<<local_index<<"]="<<use_vertex_color[local_index]<<std::endl;
           
            }

        
    }
return use_vertex_color;
}

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 std::vector< std::vector< unsigned int> > color2vertex(std::shared_ptr<dolfin::Mesh> mesh,
                                                       const MeshConnectivity &topology_N2N,
                                                       const std::map<int, std::set<unsigned int> >& shared_vertices,
                                                       unsigned int &max_num_colors)
{
int world_rank,world_size;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
std::vector< std::vector< unsigned int> > color_2_vertex;
if(world_size==1)
{
max_num_colors=2;
int num_vertices=mesh->num_vertices();
std::vector< std::vector< unsigned int> > color_2_vertex(max_num_colors);
std::vector<unsigned int> vertex_color(num_vertices);
std::vector<bool> use_vertex_color(mesh->num_vertices(),true);

color_2_vertex[0].push_back(0);
for(unsigned int ii=1;ii<num_vertices;ii++)
    color_2_vertex[1].push_back(ii);

vertex_color[0]=0;
for(unsigned int ii=1;ii<num_vertices;ii++)
    vertex_color[ii]=1;
WriteColorToFile(mesh,vertex_color,use_vertex_color);

return color_2_vertex;
}
else
{
std::vector<unsigned int> vertex_shared_global_dof;
std::vector<unsigned int> vertex_color(mesh->num_vertices(),0);

std::vector<unsigned int> used_color(1,mesh->num_vertices());
std::vector<unsigned int> keys; 
std::vector<unsigned int> shared_rank, all_shared_ranks;
std::map<unsigned int, unsigned int> shared_rank_map;
std::vector<std::set<unsigned int> > vals;

find_communicating_processes(shared_rank,all_shared_ranks, shared_vertices,shared_rank_map);
std::map<unsigned int, unsigned int> global_to_local_vertex=global2local_vertex(mesh,shared_vertices) ;
std::vector<std::vector<unsigned int> > vertex_shared_global_dof_and_color(shared_rank_map.size());



   // RECEIVE 
   for(int ii=0;ii<all_shared_ranks.size();ii++)
     {// receive from processes with a lower rank
      if(all_shared_ranks[ii]<world_rank)
         {MPI_Status status_now;
         int count_recv;
         unsigned int* buffer;
         MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status_now);
         MPI_Get_count(&status_now, MPI_UNSIGNED, &count_recv); 
         buffer=(unsigned int*)malloc(sizeof(unsigned int)*count_recv);
         MPI_Recv(buffer, count_recv, MPI_UNSIGNED, status_now.MPI_SOURCE, status_now.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);// &status_now); 
         
          for(int jj=0;jj<count_recv;jj=jj+2)
               {
               vertex_color[global_to_local_vertex.at(buffer[jj])]=buffer[jj+1];
               if(used_color.size()-1<buffer[jj])
                 while(used_color.size()<buffer[jj])
                      used_color.push_back(0);
                used_color[ buffer[jj] ] ++; 
               }
         //free (buffer);
          }
      }
 
 

 
 
        
   coloring(shared_vertices,mesh,topology_N2N,vertex_color,vertex_shared_global_dof,shared_rank_map,vertex_shared_global_dof_and_color,used_color);

    for(int ii=0;ii<all_shared_ranks.size();ii++)
      {
      if(all_shared_ranks[ii]>world_rank)
        {unsigned int* buffer;
        int tag=world_rank*world_size+all_shared_ranks[ii];
        int rank2send=shared_rank_map.at(all_shared_ranks[ii]);
        int recv_rank=all_shared_ranks[ii];
        int recv_size=vertex_shared_global_dof_and_color[rank2send].size();
        buffer=(unsigned int*)malloc(sizeof(unsigned int)*recv_size);
        for(int jj=0;jj<recv_size;jj++)
        buffer[jj]=vertex_shared_global_dof_and_color[rank2send][jj];
        MPI_Send(buffer,recv_size,MPI_UNSIGNED,recv_rank,tag,MPI_COMM_WORLD);
        }
      }


MPI_Barrier(MPI_COMM_WORLD);

unsigned int used_color_size = (*max_element(vertex_color.begin(), vertex_color.end())+1);


MPI_Allreduce(&used_color_size,&max_num_colors,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD);

std::vector< std::vector< unsigned int> > color_2_vertex(max_num_colors);



 std::vector<bool> use_vertex_color=divide_nodes_among_processes(mesh,shared_vertices,shared_rank,all_shared_ranks, shared_rank_map,global_to_local_vertex);


// --Loop on all the vertices belonging to the domain
//   Then add a vertex to the color_2_vertex map if and only if:
//   the current process has the smallest rank
//   among all the processors to which this vertex belong 
// -- If the vertex belongs only to this processor, then for sure it will be added
unsigned int num_domain_vertices=mesh->topology().ghost_offset(0);
for(unsigned int ii=0;ii<num_domain_vertices;ii++)
    {
    if(use_vertex_color[ii]==true)
    color_2_vertex[vertex_color[ii]].push_back(ii);
    }
    
WriteColorToFile(mesh,vertex_color,use_vertex_color);

return color_2_vertex;
}
}














std::vector<std::vector< int > > topologyN2N( std::shared_ptr<dolfin::Mesh> mesh, dolfin::MeshConnectivity &topology_N2F, dolfin::MeshConnectivity &topology_F2N)
{

 	std::vector<std::vector< int > > topology_N2N(mesh->num_vertices());
	unsigned int gdim = mesh->geometry().dim();  
    int world_rank,world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // loop on all the vertices belonging to the process
	for(int nn=0;nn<mesh->num_vertices();nn++)
	{
	auto actual_vertex=Vertex(*mesh, nn);
	auto point_vertex=actual_vertex.point();
	
      // loop on all the edges belonging to the vertex
	  for (MeshEntityIterator ee(  MeshEntity(*mesh, 0, nn),gdim-1 ); !ee.end(); ++ee)
	 {
		auto edge_dof=topology_N2F(nn)[ee.pos()];
		// loop on all the vertices belonging to the edge 
		for(int nn2=0;nn2<gdim;nn2++)
		   // if the node is different from the actual one, add it to its patch
		   if(nn!=topology_F2N(edge_dof)[nn2])
			topology_N2N[nn].push_back(topology_F2N(edge_dof)[nn2]);
           
	 }
	}
	return topology_N2N;
}   


std::vector<std::vector< unsigned int > > topologyN2PatchDofs( std::shared_ptr<dolfin::Mesh> mesh,std::shared_ptr<const dolfin::GenericDofMap> dofmap, dolfin::MeshConnectivity &topology_N2F)
{
	// --The domain is subdivided into world_size subdomains
	// 	 each domain has exactly num_domain_vertices vertices
	//	 but each processor will contain num_all_vertices
	// 	 because also ghost cells are present
	// --A vertex belongs to the domain if index_vertex<num_domain_vertices
	// --In topologyN2PatchDofs, all the vertices of the process are considered
	
	std::vector<std::vector< unsigned int > > Patch(mesh->num_vertices());
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
	

return Patch;
}


































PetscErrorCode smoothing_on_given_color(const std::vector< std::vector< unsigned int> > &color_2_vertex, const std::vector<std::vector< unsigned int > > &Patch,
                              Vec &corrloc, Vec &resloc, Mat* &Aloc,
                              const IS& all_local_is, const unsigned int &color)
{
PetscErrorCode ierr; 
PetscInt minusone=-1;
IS patch_is,patch_local_is;
for(std::vector<unsigned int>::const_iterator iter_c=color_2_vertex[color].begin();iter_c!=color_2_vertex[color].end();iter_c++)
   {
   
   int ii=*iter_c;   
   Mat Apatch,Apatchall; // submatrix related to the patch dofs of Aloc
   Vec respatch,corrpatch,Acorrloc; // vectors related to the patch dofs of res,x and the correction
   PetscInt sizepatch=Patch[ii].size();
   PetscInt* indices=(PetscInt*)malloc(sizepatch*sizeof(PetscInt));
   PetscInt* patch_local_indices=(PetscInt*)malloc(sizepatch*sizeof(PetscInt));
     for(int jj=0;jj<sizepatch;jj++)
        {
        indices[jj]=Patch[ii][jj];
        patch_local_indices[jj]=jj;
        }
   ierr=ISCreateGeneral(PETSC_COMM_SELF,sizepatch,indices,PETSC_COPY_VALUES,&patch_is);CHKERRQ(ierr);
   ierr=ISCreateGeneral(PETSC_COMM_SELF,sizepatch,patch_local_indices,PETSC_COPY_VALUES,&patch_local_is);CHKERRQ(ierr);
   
   ierr=MatCreateSubMatrix(*Aloc,patch_is,patch_is, MAT_INITIAL_MATRIX,&Apatch);CHKERRQ(ierr); // the patch submatrix, dim=sizepatch x sizepatch
   ierr=MatCreateSubMatrix(*Aloc,patch_is,all_local_is, MAT_INITIAL_MATRIX,&Apatchall);CHKERRQ(ierr); // the patch submatrix, dim=sizepatch x L2G_dim
   ierr=VecGetSubVector(resloc,patch_is,&respatch);CHKERRQ(ierr); // the subvector
   
   VecCreate(PETSC_COMM_SELF,&Acorrloc);
   VecSetSizes(Acorrloc,PETSC_DECIDE, sizepatch);
   VecSetType(Acorrloc,VECSEQ);
   MatMult(Apatchall,corrloc,Acorrloc);
   
   VecAXPY(respatch,minusone,Acorrloc);
   
   VecCreate(PETSC_COMM_SELF,&corrpatch);
   VecSetSizes(corrpatch,PETSC_DECIDE, sizepatch);
   VecSetType(corrpatch,VECSEQ);
   // solve patch problem
   KSP ksp_loc;
   KSPCreate(PETSC_COMM_SELF,&ksp_loc);
   KSPSetOperators(ksp_loc,Apatch,Apatch);
   KSPSolve(ksp_loc,respatch,corrpatch);
   //VecView(corrpatch, PETSC_VIEWER_STDOUT_SELF);
   KSPDestroy(&ksp_loc);
   
    
   PetscScalar* correpatchtmp;
   correpatchtmp=(PetscScalar *)malloc(sizeof(PetscScalar)*sizepatch);
   ierr=VecGetValues(corrpatch,sizepatch,patch_local_indices,correpatchtmp);CHKERRQ(ierr);
   ierr=VecSetValues(corrloc,sizepatch,indices,correpatchtmp, ADD_VALUES);CHKERRQ(ierr);
   PetscFree(indices);
   PetscFree(patch_local_indices);
   PetscFree(correpatchtmp);   
   ISDestroy(&patch_is);
   PetscFree(patch_is);
   MatDestroy(&Apatch);
   MatDestroy(&Apatchall);
   VecDestroy(&respatch);
   VecDestroy(&corrpatch);
   VecDestroy(&Acorrloc);   
}
return ierr;
}



PetscErrorCode ColoredGaussSeidel(Mat &Amat,Vec &x, Vec &bvec,Vec &Ax,Vec &res,const unsigned int &smoothing_steps,const unsigned int &max_num_colors,
                   const std::vector< std::vector< unsigned int> > &color_2_vertex, const std::vector<std::vector< unsigned int > > &Patch,
                   Vec &corrloc, Vec &resloc, Mat* &Aloc,const IS& all_local_is, const IS& all_is)
                   {
PetscErrorCode ierr; 
PetscInt minusone=-1;
PetscReal norm_res;
int world_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

VecScatter scatter;
ierr=VecScatterCreate(corrloc,NULL,x,all_is,&scatter);CHKERRQ(ierr); 

for(int ss=0;ss<smoothing_steps;ss++)
for(int cc=1;cc<max_num_colors;cc++)
{
// extract the local residual and set to zero the local correction
ierr=VecGetSubVector(res,all_is,&resloc);CHKERRQ(ierr); 
ierr=VecZeroEntries(corrloc);CHKERRQ(ierr);

//loop on internal nodes
ierr=smoothing_on_given_color(color_2_vertex,Patch,corrloc,resloc,Aloc,all_local_is,0);CHKERRQ(ierr);
//loop on boundary nodes (colored)
ierr=smoothing_on_given_color(color_2_vertex,Patch,corrloc,resloc,Aloc,all_local_is,cc);CHKERRQ(ierr);
// wait for all local corrections and update x=x+sum_i c_i
MPI_Barrier(MPI_COMM_WORLD);
ierr=VecScatterBegin(scatter,corrloc,x, ADD_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
ierr=VecScatterEnd(scatter,corrloc,x, ADD_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
// recompute the residual, Ax=A*x, res=b-A*x
ierr=MatMult(Amat,x,Ax);CHKERRQ(ierr);
ierr=VecWAXPY(res,minusone,Ax,bvec);CHKERRQ(ierr);

//compute the norm of the residual and print it
ierr=VecNorm(res, NORM_2,&norm_res);CHKERRQ(ierr);
if(world_rank==0)
{std::cout.precision(14);
std::cout<<std::scientific<<"smoothing_step=="<<ss<<", color="<<cc<<", norm_res=="<<norm_res<<std::endl;}
}

return ierr;
}








PetscErrorCode createbc(const std::vector<std::size_t> &local_to_global_map, const Mat & Amat,const Vec &bbcfoundvec,
                         IS &all_is, IS &all_local_is)
{
unsigned int cont =0;
PetscErrorCode ierr; 
PetscInt first_row; 
PetscInt last_row; 
PetscInt L2G_dim=local_to_global_map.size();
ierr=MatGetOwnershipRange(Amat,&first_row,&last_row);CHKERRQ(ierr);
PetscInt localsize=last_row-first_row;

PetscInt* one_to_L2G_dim_vec=(PetscInt *)malloc(sizeof(PetscInt)*L2G_dim);
PetscInt* one_to_localsize_vec=(PetscInt *)malloc(sizeof(PetscInt)*localsize);
PetscScalar* bc_local_dofs=(PetscScalar *)malloc(sizeof(PetscScalar)*localsize);
PetscInt * all_indices=(PetscInt *)malloc(sizeof(PetscInt)*L2G_dim);

cont=0;
for(int ii=0;ii<L2G_dim;ii++)
   {
   all_indices[ii]=local_to_global_map[ii];
   one_to_L2G_dim_vec[ii]=ii;
   if(ii<localsize)
      one_to_localsize_vec[ii]=ii+first_row;
   }



ierr=VecGetValues(bbcfoundvec,localsize,one_to_localsize_vec,bc_local_dofs);CHKERRQ(ierr);
cont=0;
for(int ii=0;ii<localsize;ii++)
{if(bc_local_dofs[ii]==1)
    {cont++;} 
}
PetscInt* bc_dofs=(PetscInt*)malloc(sizeof(PetscInt)*cont);

cont=0;
for(int ii=0;ii<localsize;ii++)
{if(bc_local_dofs[ii]==1)
    {bc_dofs[cont]=one_to_localsize_vec[ii];
     cont++;} 
}

PetscFree(bc_local_dofs);
MatZeroRows(Amat,cont,bc_dofs,1,0,0);
MatZeroRowsColumns(Amat,cont,bc_dofs,1,0,0);
ierr=ISCreateGeneral(PETSC_COMM_SELF,L2G_dim,all_indices,PETSC_COPY_VALUES,&all_is);CHKERRQ(ierr);
ierr=ISCreateGeneral(PETSC_COMM_SELF,L2G_dim,one_to_L2G_dim_vec,PETSC_COPY_VALUES,&all_local_is);CHKERRQ(ierr);

PetscFree(bc_dofs);
PetscFree(one_to_L2G_dim_vec);
PetscFree(all_indices);    


return ierr;
}


PetscErrorCode VectorToGhostedVector(Vec &vec,Vec &ghostedvec,PetscScalar* &arraylocalghosted, const IS &all_is,const std::vector<std::size_t> & local_to_global_map,
                                     const unsigned int &L2G_dim,const unsigned int &globalsize,const unsigned int &localsize )
{
PetscErrorCode ierr;
PetscInt ghostsize=L2G_dim-localsize;
arraylocalghosted=(PetscScalar *)malloc(L2G_dim*sizeof(PetscScalar)); 
PetscInt* localindex=(PetscInt *)malloc(L2G_dim*sizeof(PetscInt)); 
PetscInt* ghostindex=(PetscInt *)malloc(ghostsize*sizeof(PetscInt)); 

// create ghost solution vector 
Vec veclocalghosted;
// veclocalghosted local vector with also ghost components
ierr=VecGetSubVector(vec,all_is,&veclocalghosted);CHKERRQ(ierr);



// from the vector to the array 
for(int ii=0;ii<L2G_dim;ii++)
   localindex[ii]=ii;
VecGetValues(veclocalghosted,L2G_dim,localindex,arraylocalghosted);

for(int ii=0;ii<ghostsize;ii++)
   ghostindex[ii]=local_to_global_map[ii+localsize];

// create the global ghost vector
VecCreateGhostWithArray(PETSC_COMM_WORLD,localsize,globalsize,ghostsize,ghostindex,arraylocalghosted,&ghostedvec);
VecGhostUpdateBegin(ghostedvec,INSERT_VALUES,SCATTER_FORWARD);
VecGhostUpdateEnd(ghostedvec,INSERT_VALUES,SCATTER_FORWARD);  


PetscFree(localindex);
PetscFree(ghostindex);
VecDestroy(&veclocalghosted);

return ierr;
}



//  void coloring1(const std::map<int, std::set<unsigned int> >& shared_vertices,  std::shared_ptr<dolfin::Mesh> mesh, 
//               const std::vector<std::vector< int > >& topology_N2N,std::vector<unsigned int>& vertex_color,
//               //std::vector<unsigned int>& vertex_global_dof,std::vector<unsigned int>& vertex_shared_color,
//               std::vector<unsigned int>& vertex_shared_global_dof,const std::map<unsigned int, unsigned int>& shared_rank_map, 
//               std::vector<std::vector<unsigned int> >& vertex_shared_global_dof_and_color,std::vector<unsigned int>& used_color)
//              // std::vector<bool> &use_vertex_color)
//  {
// 
// 
// 
// unsigned int count_vertex=0;
// int world_rank,world_size;
// MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
// MPI_Comm_size(MPI_COMM_WORLD, &world_size);
// 
//    // loop on all the shared nodes
//    for (std::map< int, std::set<unsigned int> >::const_iterator shared_node=shared_vertices.begin(); shared_node!= shared_vertices.end(); ++shared_node)
//      { 
//      
//         auto shared_procs=shared_node->second;
//         std::vector<unsigned int> patch_shared_color;
//         unsigned int vertex_index=shared_node->first;
//         auto actual_vertex=Vertex(*mesh, vertex_index);
//         unsigned int global_vertex_index=actual_vertex.global_index();
//         // no previous process has colored the actual process, therefore we can add it 
//         if(vertex_color[vertex_index]==0)
//          {
//          
//          
//          if(used_color.size()==1)
//            {
//             used_color.push_back(1);
//             vertex_color[vertex_index]=1;
//            // vertex_global_dof[vertex_index]=global_vertex_index;
//             //vertex_shared_color.push_back(1);
//             vertex_shared_global_dof.push_back(global_vertex_index);
//             for ( std::set<unsigned int>::const_iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
//                  {vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(global_vertex_index);
//                   vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(vertex_color[vertex_index]);
//                  }           
//             
//             count_vertex++;
//            }
//          else
//          {
//            
//            
//       	   auto N2N=topology_N2N[vertex_index];
//      	   
//     	    // check which color I can use 
//     	    for(int ii=0;ii<N2N.size();ii++)
//    	        {
//    	          // consider all the colors >1 of the patch
//   	          if(vertex_color[N2N[ii]]>0)
//   	            patch_shared_color.push_back(vertex_color[N2N[ii]]);
//   	         }
//   	         
//   	     // sort and unique patch_shared_color
//  		 std::sort(patch_shared_color.begin(),patch_shared_color.end());
//          auto patch_shared_color_tmp = std::unique(patch_shared_color.begin(), patch_shared_color.end());
//     	 patch_shared_color.erase(patch_shared_color_tmp, patch_shared_color.end());   	
//      	 
//     	 // if the numer of colors up to now used (except the zero) is equal to the ones on the patch
//     	 // then add onother color
//     	 unsigned int used_color_size=used_color.size();
//     	 if(patch_shared_color.size()==used_color_size-1)   
//     	   { 
//     	    used_color.push_back(1);
//     	    vertex_color[vertex_index]=used_color_size;
//     	    
//     	   // vertex_global_dof[vertex_index]=global_vertex_index;
//     	   // vertex_shared_color.push_back(used_color_size);
//             vertex_shared_global_dof.push_back(global_vertex_index);
//             for ( std::set<unsigned int>::const_iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
//                  {vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(global_vertex_index);
//                   vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(vertex_color[vertex_index]);
//                  }
//             count_vertex++; 
//     	   }     
//     	 //otherwise we can opt among one of the colors already use
//     	 // we discard the ones in patch_shared_color
//     	 // and of the remaining in used_colors, we take the one which has less vertices
//     	 else
//     	   {
//     	    std::vector<unsigned int> unused_colors(used_color);
//     	    std::vector<int> unused_colors_range(unused_colors.size());
//             std::iota(unused_colors_range.begin(), unused_colors_range.end(), 0);
//             
//     	    for(int ii=0;ii<patch_shared_color.size();ii++)
//                 {unused_colors.erase( unused_colors.begin() +  patch_shared_color[ii]-ii);
//                  unused_colors_range.erase( unused_colors_range.begin() +  patch_shared_color[ii]-ii);
//                 }
//             vertex_color[vertex_index]=unused_colors_range[std::distance(unused_colors.begin(), std::min_element(unused_colors.begin(), unused_colors.end()))];
//           //  vertex_global_dof[vertex_index]=global_vertex_index;
//     	 //   vertex_shared_color.push_back(vertex_color[vertex_index]);
//             vertex_shared_global_dof.push_back(global_vertex_index);
//             count_vertex++; 
//             for ( std::set<unsigned int>::const_iterator it_procs = shared_procs.begin(); it_procs != shared_procs.end(); ++it_procs)
//                  {vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(global_vertex_index);
//                   vertex_shared_global_dof_and_color[shared_rank_map.at(*it_procs)].push_back(vertex_color[vertex_index]);
//                  }
//             }
//                                  
//           }
//           }
// 
//          patch_shared_color.clear();
//      }    
// 
// }   










// std::vector< std::vector< unsigned int> > color2vertex1(std::shared_ptr<dolfin::Mesh> mesh,
//                                                        const std::vector<std::vector< int > > &topology_N2N,
//                                                        const std::map<int, std::set<unsigned int> >& shared_vertices,
//                                                        unsigned int &max_num_colors)
// {
// int world_rank,world_size;
// MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
// MPI_Comm_size(MPI_COMM_WORLD, &world_size);
// std::vector< std::vector< unsigned int> > color_2_vertex;
// if(world_size==1)
// {
// max_num_colors=2;
// int num_vertices=mesh->num_vertices();
// std::vector< std::vector< unsigned int> > color_2_vertex(max_num_colors);
// std::vector<unsigned int> vertex_color(num_vertices);
// std::vector<bool> use_vertex_color(mesh->num_vertices(),true);
// 
// color_2_vertex[0].push_back(0);
// for(unsigned int ii=1;ii<num_vertices;ii++)
//     color_2_vertex[1].push_back(ii);
// 
// vertex_color[0]=0;
// for(unsigned int ii=1;ii<num_vertices;ii++)
//     vertex_color[ii]=1;
// WriteColorToFile(mesh,vertex_color,use_vertex_color);
// 
// return color_2_vertex;
// }
// else
// {
// std::vector<unsigned int> vertex_shared_global_dof;
// std::vector<unsigned int> vertex_color(mesh->num_vertices(),0);
// 
// std::vector<unsigned int> used_color(1,mesh->num_vertices());
// std::vector<unsigned int> keys; 
// std::vector<unsigned int> shared_rank, all_shared_ranks;
// std::map<unsigned int, unsigned int> shared_rank_map;
// std::vector<std::set<unsigned int> > vals;
// 
// find_communicating_processes(shared_rank,all_shared_ranks, shared_vertices,shared_rank_map);
// std::map<unsigned int, unsigned int> global_to_local_vertex=global2local_vertex(mesh,shared_vertices) ;
// std::vector<std::vector<unsigned int> > vertex_shared_global_dof_and_color(shared_rank_map.size());
// 
// 
// 
//    // RECEIVE 
//    for(int ii=0;ii<all_shared_ranks.size();ii++)
//      {// receive from processes with a lower rank
//       if(all_shared_ranks[ii]<world_rank)
//          {MPI_Status status_now;
//          int count_recv;
//          unsigned int* buffer;
//          MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status_now);
//          MPI_Get_count(&status_now, MPI_UNSIGNED, &count_recv); 
//          buffer=(unsigned int*)malloc(sizeof(unsigned int)*count_recv);
//          MPI_Recv(buffer, count_recv, MPI_UNSIGNED, status_now.MPI_SOURCE, status_now.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);// &status_now); 
//          
//           for(int jj=0;jj<count_recv;jj=jj+2)
//                {
//                vertex_color[global_to_local_vertex.at(buffer[jj])]=buffer[jj+1];
//                if(used_color.size()-1<buffer[jj])
//                  while(used_color.size()<buffer[jj])
//                       used_color.push_back(0);
//                 used_color[ buffer[jj] ] ++; 
//                }
//          //free (buffer);
//           }
//       }
//  
//  
// 
//  
//  
//         
//    coloring1(shared_vertices,mesh,topology_N2N,vertex_color,vertex_shared_global_dof,shared_rank_map,vertex_shared_global_dof_and_color,used_color);
// 
//     for(int ii=0;ii<all_shared_ranks.size();ii++)
//       {
//       if(all_shared_ranks[ii]>world_rank)
//         {unsigned int* buffer;
//         int tag=world_rank*world_size+all_shared_ranks[ii];
//         int rank2send=shared_rank_map.at(all_shared_ranks[ii]);
//         int recv_rank=all_shared_ranks[ii];
//         int recv_size=vertex_shared_global_dof_and_color[rank2send].size();
//         buffer=(unsigned int*)malloc(sizeof(unsigned int)*recv_size);
//         for(int jj=0;jj<recv_size;jj++)
//         buffer[jj]=vertex_shared_global_dof_and_color[rank2send][jj];
//         MPI_Send(buffer,recv_size,MPI_UNSIGNED,recv_rank,tag,MPI_COMM_WORLD);
//         }
//       }
// 
// 
// MPI_Barrier(MPI_COMM_WORLD);
// 
// unsigned int used_color_size = (*max_element(vertex_color.begin(), vertex_color.end())+1);
// 
// 
// MPI_Allreduce(&used_color_size,&max_num_colors,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD);
// 
// std::vector< std::vector< unsigned int> > color_2_vertex(max_num_colors);
// 
// 
// 
//  std::vector<bool> use_vertex_color=divide_nodes_among_processes(mesh,shared_vertices,shared_rank,all_shared_ranks, shared_rank_map,global_to_local_vertex);
// 
// 
// // --Loop on all the vertices belonging to the domain
// //   Then add a vertex to the color_2_vertex map if and only if:
// //   the current process has the smallest rank
// //   among all the processors to which this vertex belong 
// // -- If the vertex belongs only to this processor, then for sure it will be added
// unsigned int num_domain_vertices=mesh->topology().ghost_offset(0);
// for(unsigned int ii=0;ii<num_domain_vertices;ii++)
//     {
//     if(use_vertex_color[ii]==true)
//     color_2_vertex[vertex_color[ii]].push_back(ii);
//     }
//     
// WriteColorToFile(mesh,vertex_color,use_vertex_color);
// 
// return color_2_vertex;
// }
// }








































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








   
//const int* all_shared_ranks_array;
//all_shared_ranks_array = (const int*)&all_shared_ranks[0];
// MPI_Comm MPI_COMM_GHOST;
// MPI_Group MPI_GHOST_GROUP,shared_group;
// MPI_Comm_group(MPI_COMM_WORLD, &MPI_GHOST_GROUP);
// 
// MPI_Group_incl(MPI_GHOST_GROUP,all_shared_ranks.size(),all_shared_ranks_array,&shared_group);
// MPI_Comm_create(MPI_COMM_WORLD,shared_group,&MPI_COMM_GHOST);
