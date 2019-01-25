#include <dolfin.h>
#include <fstream> 

using namespace dolfin;


void WritePatchToFile(const std::string &filename,const std::vector<std::vector< unsigned int > > &Patch, const unsigned int rank)
{
std::ofstream myfile;
std::string outputname=filename+std::to_string(rank);
myfile.open (outputname);
for(int ii=0;ii<Patch.size();ii++)
   {myfile <<"Patch{";myfile<<std::to_string(ii+1);myfile<<"}=[";
    for(int jj=0;jj<Patch[ii].size();jj ++)
       {myfile << std::to_string(Patch[ii][jj]+1);myfile <<" ";}
    myfile <<"];\n";
       }
myfile <<" \n\n";
myfile.close();

}




void WriteMatrixToFile(const std::string &filename,const std::string &matrixname,const PetscInt &M,const PetscInt &N, const Mat &Amat, const int right_rank)
{

PetscErrorCode ierr; 
std::ofstream myfile;
PetscScalar *matrixtowrite;
PetscInt *M_indices;
PetscInt *N_indices;
M_indices=(PetscInt *)malloc(M*sizeof(PetscInt)); 
N_indices=(PetscInt *)malloc(N*sizeof(PetscInt)); 
matrixtowrite=(PetscScalar *)malloc(M*N*sizeof(PetscScalar)); 

for(int ii=0;ii<M;ii++)
{M_indices[ii]=ii;}//std::cout<<"M_indices="<<M_indices[ii]<<std::endl;}
for(int ii=0;ii<N;ii++)
{N_indices[ii]=ii;}//std::cout<<"N_indices="<<N_indices[ii]<<std::endl;}


Mat *Aloc;
IS is_row,is_col; // all_is=IS related to all dofs belonging to world_rank + ghost dofs and patch_is=IS related to the patch
ISCreateGeneral(PETSC_COMM_SELF,M,M_indices,PETSC_COPY_VALUES,&is_row);//CHKERRQ(ierr);
ISCreateGeneral(PETSC_COMM_SELF,N,N_indices,PETSC_COPY_VALUES,&is_col);//CHKERRQ(ierr);

MatCreateSubMatrices(Amat,1,&is_row,&is_col, MAT_INITIAL_MATRIX,&Aloc);//CHKERRQ(ierr);
//ierr=

//PetscMalloc1(n_loc,&indices); //allocate the space for vector 
//PetscMalloc1(n_loc*n_loc,&v); //allocate the space for matrix 
//PetscMalloc1(n_loc,&z); //allocate the space for vector 
MatGetValues(*Aloc,M,M_indices,N,N_indices,matrixtowrite);

if(right_rank)
{
myfile.open (filename);
myfile <<matrixname; 
myfile <<"=[";
for(int ii=0;ii<M;ii++)
   {for(int jj=0;jj<N;jj ++)
       {myfile << std::to_string(matrixtowrite[ii*N+jj]);myfile <<" ";}//std::cout<<"ii*N+jj="<<ii*N+jj<<", "<<matrixtowrite[ii*N+jj]<<std::endl;}
    myfile <<"\n";
       }
myfile <<"] \n\n";
myfile.close();
}
}


void WriteVectorToFile(const std::string &filename,const std::string &vectorname,const PetscInt &N,const Vec &vec, const int right_rank)
{

PetscErrorCode ierr; 
std::ofstream myfile;
PetscScalar *vectortowrite;
PetscInt *N_indices;
N_indices=(PetscInt *)malloc(N*sizeof(PetscInt)); 
vectortowrite=(PetscScalar *)malloc(N*sizeof(PetscScalar)); 

for(int ii=0;ii<N;ii++)
{N_indices[ii]=ii;}//std::cout<<"N_indices="<<N_indices[ii]<<std::endl;}


Vec vecloc;
IS is; // all_is=IS related to all dofs belonging to world_rank + ghost dofs and patch_is=IS related to the patch
ISCreateGeneral(PETSC_COMM_SELF,N,N_indices,PETSC_COPY_VALUES,&is);//CHKERRQ(ierr);

VecGetSubVector(vec,is,&vecloc);
//ierr=

//PetscMalloc1(n_loc,&indices); //allocate the space for vector 
//PetscMalloc1(n_loc*n_loc,&v); //allocate the space for matrix 
//PetscMalloc1(n_loc,&z); //allocate the space for vector 
VecGetValues(vecloc,N,N_indices,vectortowrite);
if(right_rank)
{
myfile.open (filename);
myfile <<vectorname; 
myfile <<"=[";
for(int ii=0;ii<N;ii++)
       {myfile << std::to_string(vectortowrite[ii]);myfile <<";\n";}
    myfile <<"\n";
myfile <<"]";
myfile.close();
}
}










 void coloring(const std::map<int, std::set<unsigned int> >& shared_vertices,  std::shared_ptr<dolfin::Mesh> mesh, 
              const std::vector<std::vector< int > >& topology_N2N,std::vector<unsigned int>& vertex_color,
              //std::vector<unsigned int>& vertex_global_dof,std::vector<unsigned int>& vertex_shared_color,
              std::vector<unsigned int>& vertex_shared_global_dof,const std::map<unsigned int, unsigned int>& shared_rank_map, 
              std::vector<std::vector<unsigned int> >& vertex_shared_global_dof_and_color,std::vector<unsigned int>& used_color)
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

   for (std::map< int, std::set<unsigned int> >::const_iterator shared_node=shared_vertices.begin(); shared_node!= shared_vertices.end(); ++shared_node)
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
//const int* all_shared_ranks_array;
//all_shared_ranks_array = (const int*)&all_shared_ranks[0];
// MPI_Comm MPI_COMM_GHOST;
// MPI_Group MPI_GHOST_GROUP,shared_group;
// MPI_Comm_group(MPI_COMM_WORLD, &MPI_GHOST_GROUP);
// 
// MPI_Group_incl(MPI_GHOST_GROUP,all_shared_ranks.size(),all_shared_ranks_array,&shared_group);
// MPI_Comm_create(MPI_COMM_WORLD,shared_group,&MPI_COMM_GHOST);

}


void WriteColorToFile(std::shared_ptr<dolfin::Mesh> mesh, const std::vector<unsigned int> &vertex_color)
{
int world_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

std::string s1,s2,s3,output_name;
s1="coloring_"; s2=std::to_string(world_rank); s3=".txt";
output_name = s1 + s2 + s3;
std::ofstream outputFile(output_name, std::ofstream::out);
for(int ii=0;ii<mesh->num_vertices();ii++)
{
    auto actual_vertex=Vertex(*mesh, ii);
 	auto point_vertex=actual_vertex.point(); 
outputFile << point_vertex[0]<<", "<<point_vertex[1]<<", "<<vertex_color[ii]<<", "<<actual_vertex.global_index()<<"\n";
    
}

}


std::vector< std::vector< unsigned int> > color2vertex(std::shared_ptr<dolfin::Mesh> mesh,
                                                       const std::vector<std::vector< int > > &topology_N2N,
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
std::vector< std::vector< unsigned int> > color_2_vertex(max_num_colors);
int num_vertices=mesh->num_vertices();
color_2_vertex[0].push_back(0);
for(unsigned int ii=1;ii<num_vertices;ii++)
    color_2_vertex[1].push_back(ii);

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
        int rank2send=shared_rank_map.at(ii);
        int recv_rank=all_shared_ranks[ii];
        int recv_size=vertex_shared_global_dof_and_color[rank2send].size();
        buffer=(unsigned int*)malloc(sizeof(unsigned int)*recv_size);
        for(int ii=0;ii<recv_size;ii++)
        buffer[ii]=vertex_shared_global_dof_and_color[rank2send][ii];
        MPI_Send(buffer,recv_size,MPI_UNSIGNED,recv_rank,tag,MPI_COMM_WORLD);
        }
      }



MPI_Barrier(MPI_COMM_WORLD);

unsigned int used_color_size = (*max_element(vertex_color.begin(), vertex_color.end())+1);


MPI_Allreduce(&used_color_size,&max_num_colors,1,MPI_UNSIGNED,MPI_MAX,MPI_COMM_WORLD);

std::vector< std::vector< unsigned int> > color_2_vertex(max_num_colors);


for(unsigned int ii=0;ii<mesh->num_vertices();ii++)
    color_2_vertex[vertex_color[ii]].push_back(ii);
    
WriteColorToFile(mesh,vertex_color);

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


std::vector<std::vector< unsigned int > > topologyN2PatchDofs( std::shared_ptr<dolfin::Mesh> mesh,std::shared_ptr<const dolfin::GenericDofMap> dofmap, dolfin::MeshConnectivity &topology_N2F)
{
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