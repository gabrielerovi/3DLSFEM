


#include <dolfin.h>
#include "MixedPoisson.h"
#include "Prova.h"
#include <fstream> 

using namespace dolfin;

  
// Source term (right-hand side)
class Source : public Expression
{
  public: 
  Source(): Expression(2) {}
  void eval(Array<double>& values, const Array<double>& x) const
  {

    values[0] = -((1.0)*( 4.0 * M_PI*M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]) - 2.0 * M_PI*M_PI *cos( M_PI * x[0]) * cos( M_PI * x[1]) ));;
    values[1] = -((1.0)*( 4.0 * M_PI*M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]) - 2.0 * M_PI*M_PI *cos( M_PI * x[0]) * cos( M_PI * x[1]) ));;
  }
};

// Boundary source for flux boundary condition
class BoundarySource : public Expression
{
public:

  BoundarySource(const Mesh& mesh) : Expression(2), mesh(mesh) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
  }

private:

  const Mesh& mesh;

};

// Sub domain for essential boundary condition
class EssentialBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS);
  }
};



int main()
{
  // Create mesh
parameters["ghost_mode"] = "shared_vertex";
parameters["reorder_vertices_gps"] = true;
parameters["reorder_cells_gps"] = true;

auto mesh = std::make_shared<UnitSquareMesh>(20,20);

const double Lambda = 1.0;
const double Mu = 1.0;
const double Beta = 1.0/(2*Mu);

auto C_equilibrium=std::make_shared<Constant>(1.0);
auto C_constitutive=std::make_shared<Constant>(1.0);
auto beta=std::make_shared<Constant>(1.0/(2*Mu));
auto alpha=std::make_shared<Constant>(((-Beta) * Lambda /(2*Lambda +2*Mu)));

  // Construct function space
  auto W = std::make_shared<MixedPoisson::FunctionSpace>(mesh);
  
  MixedPoisson::BilinearForm a(W, W);
  MixedPoisson::LinearForm L(W);






int world_rank,world_size;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
auto gdim = mesh->geometry().dim();  
std::shared_ptr<const dolfin::GenericDofMap> dofmap = W->dofmap();




//constGraph* graph;
//std::vector<ColorType> boostcolors;
//BoostGraphColoring->compute_local_vertex_coloring(graph, boostcolors);
//auto colors=mesh->color("vertex");

// if(world_rank==1)
// {
// for(int ii=0;ii<colors.size();ii++)
// {auto cell=Cell(*mesh,ii);
// std::vector<double> cell_coor;
// cell.get_vertex_coordinates(cell_coor);
// for(int jj=0;jj<0.5*cell_coor.size();jj++)
// {
// //std::cout<<"world_rank["<<world_rank<<"]: "<<"colors["<<ii<<"]: "<<colors[ii]<<" cell "<< cell_coor[2*jj]<<", "<< cell_coor[2*jj+1]<< std::endl;
// }
// 
// }
// 
// 
// }
 
// TOPOLOGICAL INFORMATION
 // node to edge 
mesh->init(0,gdim-1);
//mesh->topology().init_ghost(0,gdim-1);
auto topology_N2F=mesh->topology()(0,gdim-1);
//topology_N2F.init_ghost(0,gdim-1);
//const std::vector< std::size_t > & colors=mesh->color;
// face to node
mesh->init(gdim-1,0);
//mesh->topology().init_ghost(gdim-1,0);
MeshConnectivity topology_F2N=mesh->topology()(gdim-1,0);

// element to node
mesh->init(gdim,0);
//mesh->topology().init_ghost(gdim,0);
MeshConnectivity topology_K2N=mesh->topology()(gdim,0);

// element to face
mesh->init(gdim,gdim-1);
//mesh->topology().init_ghost(gdim,gdim-1);
MeshConnectivity topology_K2F=mesh->topology()(gdim,gdim-1);

// face to element
mesh->init(gdim-1,gdim);
//mesh->topology().init_ghost(gdim-1,gdim);
MeshConnectivity topology_F2K=mesh->topology()(gdim-1,gdim);

// node to element
mesh->init(0,gdim);
//mesh->topology().init_ghost(0,gdim);
MeshConnectivity topology_N2K=mesh->topology()(0,gdim);




std::vector<std::vector< int > > topology_N2N=topologyN2N(mesh, topology_N2F, topology_F2N);


auto coor = mesh->coordinates();


auto vertex=VertexIterator(*mesh);

//std::vector<std::list< int > > Patch(mesh->num_vertices());
//std::vector<std::list< int > > Patch(mesh->num_vertices());
//std::vector<std::vector< int > > Patch(mesh->num_vertices());



std::pair<long unsigned int, long unsigned int> ownership_range=dofmap->ownership_range();
std::vector<std::size_t> local_to_global_map;
 dofmap->tabulate_local_to_global_dofs(local_to_global_map);
 for(int ii=0;ii<local_to_global_map.size();ii++)
 {
 //std::cout<<"---xxxx------xxxx------xxxx------xxxx---- world_rank "<<world_rank<<", component="<<ii<<", local_to_global_map="<<local_to_global_map[ii]<<std::endl;
 }

//std::cout<<"--------------->world_rank "<<world_rank<<"W->dim()=="<<W->dim()<<", mesh->num_vertices() "<<mesh->num_vertices()<<" dofmap->ownership_range=["<<ownership_range.first<<", "<<ownership_range.second<<"] "<<std::endl;




auto shared_vertices=mesh->topology().shared_entities(0);
auto shared_cells = mesh->topology().shared_entities(mesh->topology().dim());
auto num_regular_vertices = mesh->topology().ghost_offset(0);

//std::cout<<" shared_cells index "<<shared_cells<<std::endl;
std::map<int, std::set<unsigned int>>::iterator it = shared_cells.begin();
std::map<int, std::set<unsigned int>>::iterator it_sharedvertex = shared_vertices.begin();


std::vector<std::vector< int > > Patch=topologyN2PatchDofs( mesh, dofmap, topology_N2F);




 
a.C_constitutive=C_constitutive;
a.C_equilibrium=C_equilibrium;
a.beta=beta;
a.alpha=alpha;

auto f = std::make_shared<Source>();
L.f = f;
L.C_equilibrium=C_equilibrium;


//PETScVector b ; 
PetscErrorCode ierr; 
PETScMatrix A;
PETScVector b;


assemble(A,a);
assemble(b,L);
Mat Amat=A.mat();
Vec bvec=b.vec();
Vec Ax,res,x;




PetscInt first_row;
PetscInt last_row;
PetscInt first_component;
PetscInt last_component;
PetscScalar *v;
PetscScalar *z;
PetscInt *indices;
PetscInt n_loc=4;
//indices=(PetscInt *)malloc(n_loc*n_loc*sizeof(PetscInt)); 
PetscMalloc1(n_loc,&indices);
PetscMalloc1(n_loc*n_loc,&v);
PetscMalloc1(n_loc,&z);


MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
//ierr=MatView(Amat, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);


MatCreateVecs(Amat,&x,&Ax);
MatMult(Amat,bvec,Ax);
ierr=MatGetOwnershipRange(Amat,&first_row,&last_row);CHKERRQ(ierr);
ierr=VecGetOwnershipRange(bvec,&first_component,&last_component);CHKERRQ(ierr);









//for(int ii=0;ii<n_loc;ii++)
//   indices[ii]=ii+first_row;
PetscInt W_dim=W->dim();
PetscInt L2G_dim=local_to_global_map.size();
PetscInt * all_indices;

all_indices=(PetscInt *)malloc(sizeof(all_indices)*L2G_dim);
for(int ii=0;ii<L2G_dim;ii++)
   {all_indices[ii]=local_to_global_map[ii];
   //std::cout<<"all_indices["<<ii<<"]="<<all_indices[ii]<<std::endl;
   }
indices[0]=0;//+first_row;
indices[1]=1;//+first_row;
indices[2]=2;//+first_row;
indices[3]=3;//+first_row;
//std::cout<<"idx "<<idx[ii]<<std::endl;
//MatGetValues(Amat,n_loc,indices,n_loc,indices,v);
//VecGetValues(bvec,n_loc,indices,z);

   
   
//Mat submat2;
//ierr=MatGetLocalSubMatrix(Amat, is, is,&submat2);CHKERRQ(ierr);
//ierr=MatRestoreLocalSubMatrix(Amat, is, is,&submat2);CHKERRQ(ierr);
//ierr=MatAssemblyEnd(submat2,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

Mat* Aloc; // the diagonal submatrix related to all dofs belonging to world_rank + ghost dofs
Vec xloc,bloc;  // the subvector related to all dofs belonging to world_rank + ghost dofs
VecCreate(PETSC_COMM_SELF,&xloc);
VecCreateSeq(PETSC_COMM_SELF ,L2G_dim,&xloc);
VecZeroEntries(xloc);

IS all_is, patch_is; // IS related to all dofs belonging to world_rank + ghost dofs and IS related to the patch
ierr=ISCreateGeneral(PETSC_COMM_SELF,n_loc,indices,PETSC_COPY_VALUES,&patch_is);CHKERRQ(ierr);
ierr=ISCreateGeneral(PETSC_COMM_SELF,L2G_dim,all_indices,PETSC_COPY_VALUES,&all_is);CHKERRQ(ierr);
//ierr=ISView(is,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);


ierr=VecGetSubVector(bvec,all_is,&bloc);CHKERRQ(ierr); // the subvector
//ierr=VecView(bloc, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector
//ierr=VecRestoreSubVector(bvec,is,&bloc);CHKERRQ(ierr);

ierr=MatCreateSubMatrices(Amat,1,&all_is,&all_is, MAT_INITIAL_MATRIX,&Aloc);CHKERRQ(ierr); // the diagonal submatrix
//ierr=MatView(*Aloc, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the diagonal submatrix



Mat Apatch;
Vec respatch,xpatch,corrpatch;
PetscInt patchsize;
ierr=MatCreateSubMatrix(*Aloc,patch_is,patch_is, MAT_INITIAL_MATRIX,&Apatch);CHKERRQ(ierr); // the patch submatrix
//ierr=MatView(Apatch, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the patch submatrix
ierr=VecGetSubVector(bloc,patch_is,&respatch);CHKERRQ(ierr); // the subvector
//ierr=VecView(respatch, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the patch residual
 
 
VecGetSize(respatch,&patchsize);
VecCreate(PETSC_COMM_SELF,&corrpatch);
VecSetSizes(corrpatch,PETSC_DECIDE, patchsize);
VecSetType(corrpatch,VECSEQ);

KSP ksp_loc;
KSPCreate(PETSC_COMM_SELF,&ksp_loc);
KSPSetOperators(ksp_loc,Apatch,Apatch);
KSPSolve(ksp_loc,respatch,corrpatch);
//VecView(corrpatch, PETSC_VIEWER_STDOUT_SELF);
KSPDestroy(&ksp_loc);


PetscScalar* correpatchtmp;
correpatchtmp=(PetscScalar *)malloc(sizeof(PetscScalar)*patchsize);
ierr=VecGetValues(corrpatch,patchsize,indices,correpatchtmp);CHKERRQ(ierr);
ierr=VecSetValues(xloc,patchsize,indices,correpatchtmp, ADD_VALUES);CHKERRQ(ierr);
//ierr=VecView(xloc, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the patch residual


PetscFree(indices);
PetscFree(all_indices);
ISDestroy(&patch_is);
PetscFree(patch_is);
ISDestroy(&all_is);
PetscFree(all_is);

////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////           FIND COMMUNICATING PROCESSES           /////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<unsigned int> keys; 
std::vector<unsigned int> shared_rank, all_shared_ranks;
std::map<unsigned int, unsigned int> shared_rank_map;
std::vector<std::set<unsigned int> > vals;
const int* all_shared_ranks_array;
unsigned int cont=0;
     
    for (std::map< int, std::set<unsigned int> >::iterator it=shared_vertices.begin(); it!= shared_vertices.end(); ++it)
      {keys.push_back(it->first);
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
all_shared_ranks_array = (const int*)&all_shared_ranks[0];

MPI_Comm MPI_COMM_GHOST;
MPI_Group MPI_GHOST_GROUP,shared_group;
MPI_Comm_group(MPI_COMM_WORLD, &MPI_GHOST_GROUP);

MPI_Group_incl(MPI_GHOST_GROUP,all_shared_ranks.size(),all_shared_ranks_array,&shared_group);
MPI_Comm_create(MPI_COMM_WORLD,shared_group,&MPI_COMM_GHOST);


for(int ii=0;ii<all_shared_ranks.size();ii++)
   {
   //std::cout<<"world_rank=="<<world_rank<<", shared_rank : "<<shared_rank[ii]<<std::endl;
   }

//// Create Map Global to Local for vertices. It is important because the global indices are fixed and we use them for coloring
std::map<unsigned int, unsigned int> global_to_local_vertex;   
 for (std::map< int, std::set<unsigned int> >::iterator shared_node=shared_vertices.begin(); shared_node!= shared_vertices.end(); ++shared_node)
 {
   auto vertex_index=shared_node->first;
   auto actual_vertex=Vertex(*mesh, vertex_index);
   global_to_local_vertex[actual_vertex.global_index()] = vertex_index;
    
 }

std::map<unsigned int,unsigned int>::iterator it_global_to_local_vertex = global_to_local_vertex.begin();
while(it_global_to_local_vertex != global_to_local_vertex.end())
{
//std::cout<<"world_rank: "<<world_rank<<"it_global_to_local_vertex: "<<it_global_to_local_vertex->first<<" :: "<<global_to_local_vertex[it_global_to_local_vertex->first]<<std::endl;
it_global_to_local_vertex++;
}

 
 
 
 
 
 
 

// we initialize with zero the vector, since 0 is the color related to non-shared (internal) nodes
unsigned int num_shared_vertices= std::distance(shared_vertices.begin(),shared_vertices.end());

     
     
  
std::vector<unsigned int> vertex_shared_color;
std::vector<unsigned int> vertex_shared_global_dof;
std::vector<unsigned int> vertex_color(mesh->num_vertices(),0);
std::vector<unsigned int> vertex_global_dof(mesh->num_vertices(),0);

for(int ii=0;ii<shared_rank.size();ii++)
   shared_rank_map[shared_rank[ii]]=ii;

std::vector<std::vector<unsigned int> > vertex_shared_global_dof_and_color(shared_rank_map.size());
//coloring(shared_vertices,mesh,topology_N2N,vertex_color, vertex_global_dof,vertex_shared_color,vertex_shared_global_dof,shared_rank_map,vertex_shared_global_dof_and_color);
// for(int jj=0;jj<shared_rank.size();jj++)
// for(int ii=0;ii<vertex_shared_global_dof_and_color[jj].size();ii=ii+2)
//       std::cout<<"SHARED world_rank: "<<world_rank<<" related to "<<shared_rank[jj]<<"-- "<<vertex_shared_global_dof_and_color[jj][ii]<<", " <<vertex_shared_global_dof_and_color[jj][ii+1]<<std::endl;

     

std::vector<MPI_Request> request;
std::vector<MPI_Status> status;
std::vector<int> ricevimi(world_size,0);
std::vector<int> mandami(world_size,0);
 
 
 
 
 
  
std::string s1,s2,s3,output_name;
s1="sending_"; s2=std::to_string(world_rank); s3=".txt";
output_name = s1 + s2 + s3;
std::ofstream outputFileSend(output_name, std::ofstream::out);
































MPI_Barrier(MPI_COMM_WORLD);
std::vector<unsigned int> used_color(1,mesh->num_vertices());






   
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
         std::cout<<"world_rank=="<<world_rank<<" PRE RECV, source=="<<status_now.MPI_SOURCE<<" , tag=="<<status_now.MPI_TAG<<std::endl;
         MPI_Recv(buffer, count_recv, MPI_UNSIGNED, status_now.MPI_SOURCE, status_now.MPI_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);// &status_now); 
         std::cout<<"world_rank=="<<world_rank<<" POST RECV, source=="<<status_now.MPI_SOURCE<<" , tag=="<<status_now.MPI_TAG<<std::endl;
         
          for(int jj=0;jj<count_recv;jj=jj+2)
               {vertex_color[global_to_local_vertex[buffer[jj]]]=buffer[jj+1];
               //std::cout<<" RICEVO world_rank: "<<world_rank<<" count_recv   "<<count_recv<<", "<<jj<<", ("<<global_to_local_vertex[buffer[jj]]<<", "<<buffer[jj]<<", "<<buffer[jj+1]<<")"<<std::endl;
               if(used_color.size()-1<buffer[jj])
                 while(used_color.size()<buffer[jj])
                      used_color.push_back(0);
                used_color[ buffer[jj] ] ++; 
               }
         //free (buffer);
          }
      }
        
   coloring(shared_vertices,mesh,topology_N2N,vertex_color, vertex_global_dof,vertex_shared_color,vertex_shared_global_dof,shared_rank_map,vertex_shared_global_dof_and_color,used_color);
std::cout<<"world_rank=="<<world_rank<<" post coloring out "<<std::endl;

    for(int ii=0;ii<all_shared_ranks.size();ii++)
      {
      std::cout<<"world_rank=="<<world_rank<<" ii "<<ii<<"/"<<all_shared_ranks.size()<<std::endl;
      if(all_shared_ranks[ii]>world_rank)
        {unsigned int* buffer;
        int tag=world_rank*world_size+all_shared_ranks[ii];
        int rank2send=shared_rank_map[ii];
        int recv_rank=all_shared_ranks[ii];
        int recv_size=vertex_shared_global_dof_and_color[rank2send].size();
        //=&vertex_shared_global_dof_and_color[rank2send][0];
        buffer=(unsigned int*)malloc(sizeof(unsigned int)*recv_size);
        for(int ii=0;ii<recv_size;ii++)
	      buffer[ii]=vertex_shared_global_dof_and_color[rank2send][ii];
	      
	    std::cout<<"world_rank=="<<world_rank<<" PRE SEND, source=="<<recv_rank<<" , tag=="<<tag<<" all_shared_ranks"<<all_shared_ranks[ii]<<std::endl;
        MPI_Send(buffer,recv_size,MPI_UNSIGNED,recv_rank,tag,MPI_COMM_WORLD);
        std::cout<<"world_rank=="<<world_rank<<" POST SEND, source=="<<recv_rank<<" , tag=="<<" all_shared_ranks"<<all_shared_ranks[ii]<<tag<<std::endl;

    //     for(int jj=0;jj<recv_size;jj++)
//         std::cout<<"MANDA world_rank: "<<world_rank<<" to "<<all_shared_ranks[ii]<<" tag=="<<tag<<". dimensione  "<<vertex_shared_global_dof_and_color[rank2send].size()<<", "<<vertex_shared_global_dof_and_color[rank2send][jj]<<", buffer "<<buffer[jj]<<std::endl;
        }
        }






MPI_Barrier(MPI_COMM_WORLD);
// 
// int rank2send=shared_rank_map[1];
// int recv_rank=all_shared_ranks[1];
// int MAX_NUMBERS = vertex_shared_global_dof_and_color[rank2send].size();
// unsigned int* numbers;
// int number_amount;
// 
// if (world_rank == 0) {
//     // Pick a random amount of integers to send to process one
//     srand(time(NULL));
//     number_amount = MAX_NUMBERS ;
//     numbers=(unsigned int *)malloc(MAX_NUMBERS *sizeof(unsigned int));
// 
//     for(int ii=0;ii<MAX_NUMBERS;ii++)
//      numbers[ii]=vertex_shared_global_dof_and_color[rank2send][ii];
//     MPI_Send(numbers, number_amount, MPI_UNSIGNED, 1, 0, MPI_COMM_WORLD);
// } else if (world_rank == 1) {
//     MPI_Status status;
//     
//     MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//     MPI_Get_count(&status, MPI_UNSIGNED, &number_amount); 
//     numbers=(unsigned int *)malloc(number_amount *sizeof(unsigned int));
//     MPI_Recv(numbers, number_amount, MPI_UNSIGNED, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//     for(int ii=0;ii<number_amount;ii++)
//      std::cout<<" recv num "<<numbers[ii]<<std::endl;
//            
// }
















s1="coloring_"; s2=std::to_string(world_rank); s3=".txt";
output_name = s1 + s2 + s3;
std::ofstream outputFile(output_name, std::ofstream::out);
for(int ii=0;ii<mesh->num_vertices();ii++)
{
    auto actual_vertex=Vertex(*mesh, ii);
 	auto point_vertex=actual_vertex.point(); 
outputFile << point_vertex[0]<<", "<<point_vertex[1]<<", "<<vertex_color[ii]<<", "<<actual_vertex.global_index()<<"\n";
    
}
for(int jj=0;jj<world_size;jj++)
     {
      MPI_Barrier(MPI_COMM_WORLD);
     if(world_rank==jj)
          {
          for(int ii=0;ii<keys.size();ii++)
			{
 			// std::cout<<"KEYS NUMBER: "<<ii<<std::endl;
 			 auto vertex_index=keys[ii];
 			 auto actual_vertex=Vertex(*mesh, vertex_index);
 			 auto point_vertex=actual_vertex.point();

 	       //std::cout<<"world_rank: "<<world_rank<<"keys "<<keys[ii]<<"  point=(" <<point_vertex[0]<<", "<<point_vertex[1]<<") "<<std::endl;
        	for ( std::set<unsigned int>::iterator setit = vals[ii].begin(); setit != vals[ii].end(); ++setit)
        	     {
        	     }
             //std::cout<<"--- "<<*setit<<std::endl; 
             } 
             }
	  }
	 


  for(int jj=0;jj<world_size;jj++)
     {
     MPI_Barrier(MPI_COMM_WORLD);
     if(world_rank==jj)
while(it != shared_cells.end())
{
auto nodes=topology_K2N(it->first);
//std::cout<<"world_rank: "<<world_rank<<", cell "<<it->first<<", nodes "<< nodes[0]<<", "<<nodes[1]<<", "<<nodes[2]<<std::endl;
for(int ii=0;ii<3;ii++)
{
auto actual_vertex=Vertex(*mesh, nodes[ii]);
auto point_vertex=actual_vertex.point();
//std::cout<<point_vertex[0]<<", "<<point_vertex[1]<<std::endl;
}
it++;
}
}


// MatView(submat, PETSC_VIEWER_STDOUT_WORLD);

//assemble(b,L);


// Define boundary condition
auto G = std::make_shared<BoundarySource>(*mesh);
auto boundary = std::make_shared<EssentialBoundary>();
DirichletBC bc(W->sub(1), G, boundary);


  // Compute solution
  Function w(W);
  solve(a == L, w, bc);

// Now, the separate components ``sigma`` and ``u`` of the solution can
// be extracted by taking components. These can then be written to XDMF
// files.
// 
// .. code-block:: cpp

  // Extract sub functions (function views)
   Function& sigma = w[0];
   Function& u = w[1];
 
   File file1("mixedpoisson_u.pvd");
   file1 << u;
   File file2("mixedpoisson_sigma.pvd");
   file2 << sigma;
  
//   XDMFFile("sigma.xdmf").write(sigma);
//   XDMFFile("u.xdmf").write(u);

  return 0;
}

















//MeshTopology topology=MeshTopology();
// loop on all the vertices of the local 
//for (; !vertex.end(); ++vertex)
//{ 
// auto prova=vertex->point();
// std::vector<long unsigned int> vertex_index(1);
// vertex_index[0]=vertex->index();
// 
// const unsigned int* vertex_to_faces=topology_N2F(vertex_index[0]);
// auto tmp=dofmap->entity_dofs(*mesh, 0,vertex_index);
// 
// //std::cout<<"world_rank: "<<world_rank<<", vertex: "<<vertex->index()<<"ok "<<tmp.size()<<std::endl;
// //std::cout<<"world_rank: "<<world_rank<<", vertex: "<<prova[0]<<", "<<prova[1]<<std::endl;
// //for(int ii=0;ii<sizeof(ww)/sizeof(ww[0]);ii++)
//  //std::cout<<"world_rank: "<<world_rank<<", vertex: "<<vertex->index()<<" ii= "<<ii<< ", dofs "<<tmp[ii]<< std::endl;   
// 
// std::vector<long unsigned int> N2F(vertex_to_faces,vertex_to_faces+sizeof(vertex_to_faces)/sizeof(vertex_to_faces[0]));
// //test(prova2)
// //std::copy(vertex_to_faces,vertex_to_faces+sizeof(vertex_to_faces)/sizeof(vertex_to_faces[0]),prova2.begin());
// auto tmp_face=dofmap->entity_dofs(*mesh, gdim-1,N2F);
// auto tmp_node=dofmap->entity_dofs(*mesh, 0,vertex_index);
// 
// if(world_rank==1)
// std::cout<<"world_rank: "<<world_rank<<"not shared vertex_index " <<coor[gdim*vertex_index[0]]<<", "<<coor[1+gdim*vertex_index[0]]<<", "<<vertex_index[0]<<std::endl;
// 
// for(int ii=0;ii<N2F.size();ii++)
// {
// //std::cout<<"world_rank: "<<world_rank<<"not shared vertex_index " <<coor[gdim*vertex_index[0]]<<", "<<coor[1+gdim*vertex_index[0]]<<", "<<vertex_index[0]<<" N2F " <<N2F[ii]<<std::endl;
// }
// 
// for(int ii=0; ii<sizeof(vertex_to_faces)/sizeof(vertex_to_faces[0]);ii++)
// {
// //std::cout<<"world_rank: "<<world_rank<<", vertex_to_faces: "<<vertex_to_faces[ii]<<std::endl;
// }
// 
// for(int ii=0; ii<tmp_face.size();ii++)
// {
// Patch[vertex_index[0]].push_back(tmp_face[ii]);
// }
// 
// 
// for(int ii=0; ii<gdim;ii++)
// {
// Patch[vertex_index[0]].push_back(tmp_node[ii]);
// }
// 
// // if(world_rank==2)
// // {
// // for(int ii=0; ii<(tmp_face.size()+gdim);ii++)
// // {
// // std::cout<<"world_rank: "<<world_rank<<", vertex_index: "<<vertex_index[0]<<", patch: "<<Patch[vertex_index[0]][ii]<<std::endl;
// // }
// // }
// 
//}

// for(int ii=0;ii<5;ii++)
//  for(int jj=0;jj<sizeof(topology_F2N(ii))/sizeof(topology_F2N(ii)[0]);jj++)
// {}//std::cout<<"topology_F2N " <<ii<< ", "<<jj<<", "<<topology_F2N(ii)[jj]<<std::endl;
// 
// for(int ii=0;ii<4;ii++)
// {
// int cont=0;
// //auto ent=MeshEntity(*mesh, 0, ii);
// //MeshEntityIterator ee( MeshEntity(*mesh, 0, ii),gdim-1 );
// 
// for (MeshEntityIterator ee( MeshEntity(*mesh, 0, ii),gdim-1 ); !ee.end(); ++ee)
// {
//   
//   //std::cout<<"topology_N2F" <<ii<< ", "<<cont<<", "<<topology_N2F(ii)[cont]<<std::endl;
//   cont++;
// }
// }
// 
// 
// if(world_rank==0)
// {
// 
// while(it != shared_cells.end())
// {
// 
// auto dof_c = dofmap->cell_dofs(it->first);
// 
// for(int ii=0;ii<dof_c.size();ii++)
// {
// //printf("shared dof_c= %d\n ", dof_c[ii]);
// }
//  for(int ii=0;ii<3;ii++)
//  {
// 
// int p=topology_K2N(it->first)[ii];
// std::vector<long unsigned int> vertex_index(1);
// vertex_index[0]=p;
// const unsigned int* vertex_to_faces=topology_N2F(p);
// std::vector<long unsigned int> N2F(vertex_to_faces,vertex_to_faces+sizeof(vertex_to_faces)/sizeof(vertex_to_faces[0]));
// auto tmp_face=dofmap->entity_dofs(*mesh, gdim-1,N2F);
// auto tmp_node=dofmap->entity_dofs(*mesh, 0,vertex_index);
// //std::cout<<"SHARED world_rank: "<<world_rank<<"topology_K2N: "<<it->first<<", "<<p<<std::endl;
// //std::cout<<" SHARED world_rank: "<<world_rank<<"cell: "<<it->first<<", "<<p<<std::endl;
// 
// 
// for(int ii=0;ii<N2F.size();ii++)
// {
// //std::cout<<" vertex_index " <<p<<" N2F " <<N2F[ii]<<std::endl;
// }
// for(int ii=0; ii<tmp_face.size();ii++)
// {
// //std::cout<<" vertex_index " <<p<<" tmp_face " <<tmp_face[ii]<<std::endl;
// //Patch[vertex_index[0]].push_back(tmp_face[ii]);
// }
// 
// for(int ii=0; ii<gdim;ii++)
// {
// //std::cout<<" tmp_node" <<tmp_node[ii]<<std::endl;
// //Patch[vertex_index[0]].push_back(tmp_node[ii]);
// }
// 
// 
// 
// 
// 
// }
// //std::cout<<std::endl<<"world_rank cell: "<<world_rank<<"it first"<<it->first<<std::endl;
// //for (std::set<unsigned int>::iterator itvalue=it->second.begin(); itvalue!=it->second.end(); ++itvalue)
// //	    std::cout <<' '<< *itvalue;
// 
// it++;
// }
// 
// 
// auto cell=CellIterator(*mesh);
// // loop on all the vertices of the local 
// for (; !cell.end(); ++cell)
// {
// 
// auto dof_c = dofmap->cell_dofs(cell->index());
// 
// for(int ii=0;ii<dof_c.size();ii++)
// {
// //printf("dof_c= %d\n ", dof_c[ii]);
// }
// 
// 
// for(int ii=0;ii<3;ii++)
// {
// int p=topology_K2N(cell->index())[ii];
// std::vector<long unsigned int> vertex_index(1);
// vertex_index[0]=p;
// //std::cout<<"NOT SHARED world_rank: "<<world_rank<<"cell: "<<cell->index()<<", "<<p<<std::endl;
// const unsigned int* vertex_to_faces=topology_N2F(p);
// std::vector<long unsigned int> N2F(vertex_to_faces,vertex_to_faces+sizeof(vertex_to_faces)/sizeof(vertex_to_faces[0]));
// auto tmp_face=dofmap->entity_dofs(*mesh, gdim-1,N2F);
// auto tmp_node=dofmap->entity_dofs(*mesh, 0,vertex_index);
// 
// //for(int ii=0; ii<sizeof(vertex_to_faces)/sizeof(vertex_to_faces[0]);ii++)
// for(int ii=0; ii<3;ii++)
// {
// //std::cout<<sizeof(topology_N2F)<<", "<<sizeof(topology_N2F(0))<<", "<<sizeof(topology_N2F(1))<<", "<<sizeof(topology_N2F(2))<<", "<<sizeof(topology_N2F(3))<<", "<<sizeof(topology_N2F(0)[0])<<", "<<p<<" vertex_to_faces " <<vertex_to_faces[ii]<<std::endl;
// //Patch[vertex_index[0]].push_back(tmp_face[ii]);
// }
// 
// 
// for(int ii=0; ii<gdim;ii++)
// {
// //std::cout<<" tmp_node" <<tmp_node[ii]<<std::endl;
// //Patch[vertex_index[0]].push_back(tmp_node[ii]);
// }
// }
// 
// }
// 
// }
//std::cout<<" world_rank: "<<world_rank<<" rbegin()->first; " << shared_cells.rbegin()->first<<std::endl;

//for(int ii=0;ii<coor.size();ii++)
//   std::cout<<"world_rank "<<world_rank<<", ii= "<<ii<<", "<<coor[ii]<<std::endl;

//if(world_rank==0)

//it=shared_vertices.begin();
// for (std::map<unsigned int, std::set<int> >::iterator itvalue=shared_vertices.begin(); itvalue!=shared_vertices.end(); ++itvalue)
// {
// }
 // while(it != shared_vertices.end())
//  {
// //std::cout<<std::endl<<"world_rank: "<<world_rank<<"vertex"<<it->first<<std::endl;
// for (std::set<unsigned int>::iterator itvalue=it->second.begin(); itvalue!=it->second.end(); ++itvalue)
//      {}//std::cout <<' '<<*itvalue<<' '<< coor[gdim*(it->first)]<<", "<< coor[gdim*(it->first)]<<std::endl; 
//  it++;
//  }




//auto coordinates=mesh->coordinates();
//     std::vector<unsigned int> keys;
//     std::vector<std::set<unsigned int> > vals;
//     for (std::map< int, std::set<unsigned int> >::iterator it=shared_vertices.begin(); it!= shared_vertices.end(); ++it)
//        {keys.push_back(it->first);
//        vals.push_back(it->second);
//        }
// 
// for(int ii=0;ii<keys.size();ii++)
//  {
//  //std::cout<<"world_rank: "<<world_rank<<"keys &val  "<<keys[ii]<<std::endl;
// for ( std::set<unsigned int>::iterator setit = vals[ii].begin(); setit != vals[ii].end(); ++setit)
// {
// //std::cout<<" val  "<<*setit<<std::endl;
// //std::cout<<coordinates[gdim*ii]<<", "<<coordinates[gdim*ii+1]<<std::endl;
// }
//  
//  }
// std::cout<<"oltre "<<coordinates[gdim*2]<<", "<<coordinates[gdim*2+1]<<std::endl;
//  std::cout<<"oltre "<<coordinates[gdim*6]<<", "<<coordinates[gdim*6+1]<<std::endl;
//std::cout<<"oltre 3 "<<coordinates.size()<<std::endl;
//auto keys=get_keys(shared_vertices);
//auto values=get_vals(shared_vertices);



//auto iteratore=num_regular_vertices.begin();

//std::cout<<std::endl<<"world_rank: "<<world_rank<<" num_regular_vertices "<<num_regular_vertices<<std::endl;

