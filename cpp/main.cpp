


#include <dolfin.h>
#include "MixedPoisson.h"
#include "Prova.h"
#include <fstream> 

using namespace dolfin;

// Then follows the definition of the coefficient functions (for
// :math:`f` and :math:`G`), which are derived from the DOLFIN
// :cpp:class:`Expression` class.
// 
// .. code-block:: cpp



//   std::vector<unsigned int> get_keys(std::map<  int,std::set<unsigned int > > &themap )
//   {
//     std::vector<unsigned int> keys;
//     for (std::map<unsigned int, std::set<unsigned int> >::iterator it=themap.begin(); it!= themap.end(); it++)
//        keys.push_back(it->first);
//     return keys;
//   }
//   
//   
//     std::vector<unsigned int> get_vals(std::map<unsigned int,std::set<unsigned int > >& themap , unsigned int i)
//   {
//     std::vector<unsigned int> vals;
//     for (std::set<unsigned int>::iterator it=themap[i].begin(); it!= themap[i].end(); it++)
//       vals.push_back(*it);
//     return vals;
//   } 
//   
  
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


// Then follows the definition of the essential boundary part of the
// boundary of the domain, which is derived from the
// :cpp:class:`SubDomain` class.
// 
// .. code-block:: cpp

// Sub domain for essential boundary condition
class EssentialBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS);
  }
};

// Inside the ``main()`` function we first create the ``mesh`` and then
// we define the (mixed) function space for the variational
// formulation. We also define the bilinear form ``a`` and linear form
// ``L`` relative to this function space.
// 
// .. code-block:: cpp

int main()
{
  // Create mesh
parameters["ghost_mode"] = "shared_vertex";
parameters["reorder_vertices_gps"] = true;
parameters["reorder_cells_gps"] = true;

auto mesh = std::make_shared<UnitSquareMesh>(9,9);

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
auto dofmap = W->dofmap();




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
auto topology_F2N=mesh->topology()(gdim-1,0);

// element to node
mesh->init(gdim,0);
//mesh->topology().init_ghost(gdim,0);
auto topology_K2N=mesh->topology()(gdim,0);

// element to face
mesh->init(gdim,gdim-1);
//mesh->topology().init_ghost(gdim,gdim-1);
auto topology_K2F=mesh->topology()(gdim,gdim-1);

// face to element
mesh->init(gdim-1,gdim);
//mesh->topology().init_ghost(gdim-1,gdim);
auto topology_F2K=mesh->topology()(gdim-1,gdim);

// node to element
mesh->init(0,gdim);
//mesh->topology().init_ghost(0,gdim);
auto topology_N2K=mesh->topology()(0,gdim);





std::vector<std::vector< int > > topology_N2N(mesh->num_vertices());

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







auto coor = mesh->coordinates();


auto vertex=VertexIterator(*mesh);

//std::vector<std::list< int > > Patch(mesh->num_vertices());
//std::vector<std::list< int > > Patch(mesh->num_vertices());
std::vector<std::vector< int > > Patch(mesh->num_vertices());
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



// loop on all the vertices (not ghost, but the patch will consider also the ghost faces)
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
while(it_sharedvertex != shared_vertices.end())
{
auto vertex_index=it_sharedvertex->first;
auto actual_vertex=Vertex(*mesh, vertex_index);
auto point_vertex=actual_vertex.point();

for (MeshEntityIterator ee(  MeshEntity(*mesh, 0, it_sharedvertex->first),gdim-1 ); !ee.end(); ++ee)
{
  auto edge_num=topology_N2F(it_sharedvertex->first)[ee.pos()];
  auto edge=Edge(*mesh,edge_num);
  
  auto edge_dof=topology_N2F(vertex_index)[ee.pos()];
  auto actual_edge=Edge(*mesh,edge_dof);
  auto point_edge=actual_edge.midpoint();
 // std::cout<<"SHARED world_rank: "<<world_rank<<", topology_N2F: " <<it_sharedvertex->first<<" coord "<<point_vertex[0]<<", "<<point_vertex[1]<< ", edge: "<<edge_dof<<" coord "<< point_edge[0]<<", "<<point_edge[1] <<std::endl;
}
it_sharedvertex++;
}


// loop on the vertices
for(int nn=0;nn<mesh->num_vertices() ;nn++)
    {// loop on the dofs of the given vertex
    auto actual_vertex=Vertex(*mesh, nn);
    auto point_vertex=actual_vertex.point();
    for(int ii=0;ii<Patch[nn].size();ii++)
        {
        // std::cout<<"NOT SHARED["<<world_rank<<"], vertex="<<nn<<",with coord["<< point_vertex[0]<<", "<<point_vertex[1]<<"], and dof["<<ii<<"]="<<Patch[nn][ii]<<std::endl;
        }
    }

















 
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

 
 
 
 
 
 
 

std::vector<unsigned int> used_color(1,mesh->num_vertices());
// we initialize with zero the vector, since 0 is the color related to non-shared (internal) nodes
unsigned int num_shared_vertices= std::distance(shared_vertices.begin(),shared_vertices.end());


std::vector<unsigned int> vertex_color(mesh->num_vertices(),0);//num_shared_vertices,0 );
std::vector<unsigned int> vertex_global_dof(mesh->num_vertices(),0);//num_shared_vertices,0 ); 


// for(int ii=0;ii<vertex_shared_color.size();ii++)
//     {
//     auto actual_vertex=Vertex(*mesh, ii);
//  	auto point_vertex=actual_vertex.point();    
//    // std::cout<<"//////////////world_rank: "<<world_rank<<"  point=(" <<point_vertex[0]<<", "<<point_vertex[1]<<") "<<" vertex_shared_color["<<ii<<"]=="<<vertex_shared_color[ii]<<std::endl;
//      }
     
     
     
     
     
     
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
////////////////////////////////////////             COLORING                  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     



// loop on all the processors of the world_rank processor, communicating through ghost cells
for(std::vector<unsigned int>::iterator rank_shared_on_this_proc=all_shared_ranks.begin();rank_shared_on_this_proc!=all_shared_ranks.end();rank_shared_on_this_proc++)
{
std::vector<unsigned int> vertex_shared_color;
std::vector<unsigned int> vertex_shared_global_dof;//num_shared_vertices,0 ); 
unsigned int count_vertex=0;
// consider the rank=rank_now
unsigned int rank_now=*rank_shared_on_this_proc;
// then only the rank_now processor will color the vertices, while all the others will wait (MPI_Barrier) for it to finish
// then rank_now will send to all the shared processor the information regarding the colored vertices
// so that each of them can update their color
if(rank_now==world_rank)
  {
     // loop on all the shared (ghost) vertices and use the coloring algorithm    
     for (std::map< int, std::set<unsigned int> >::iterator shared_node=shared_vertices.begin(); shared_node!= shared_vertices.end(); ++shared_node)
        {
        std::vector<unsigned int> patch_shared_color;
        unsigned int vertex_index=shared_node->first;
        auto actual_vertex=Vertex(*mesh, vertex_index);
        unsigned int global_vertex_index=actual_vertex.global_index();
         // for the first node, use the color=1
         if(used_color.size()==1)
           {
            used_color.push_back(1);
            vertex_color[vertex_index]=1;
            vertex_global_dof[vertex_index]=global_vertex_index;
            vertex_shared_color.push_back(1);
            vertex_shared_global_dof.push_back(global_vertex_index);
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
            }
                                 
         }
        patch_shared_color.clear();
        
        }
      }
      
  //  std::cout<<"PRE world_rank: "<<world_rank<<" vertex_shared_color.size() "<<vertex_shared_color.size()<<std::endl;
    for(int ii=0;  ii<vertex_shared_color.size();ii++)
       {}//std::cout<<"PRE world_rank: "<<world_rank<<" vertex_shared_color "<<vertex_shared_color[ii]<<std::endl;
 //    MPI_Bcast(&vertex_shared_color[0],count_vertex,MPI_UNSIGNED,rank_now,MPI_COMM_GHOST);
//     MPI_Bcast(&vertex_shared_global_dof[0],count_vertex,MPI_UNSIGNED,rank_now,MPI_COMM_GHOST);
//     MPI_Bcast(&count_vertex,1,MPI_UNSIGNED,rank_now,MPI_COMM_GHOST);
//     
//    MPI_Barrier(MPI_GHOST_GROUP);     
//    std::cout<<"world_rank=="<<world_rank<<" rank_now=="<< rank_now<<" count_vertex=="<<count_vertex<<std::endl;
    //std::cout<<"world_rank: "<<world_rank<<" vertex_shared_color "<<vertex_shared_color[0]<<std::endl;
    if(world_rank!=rank_now)
    {
      for(int ii=0;ii<count_vertex;ii++)
      {
//        // check whether the global dof from rank_now also belongs to world_rank
//        // in this case, update the color
       // {std::cout<<"world_rank: "<<world_rank<<" vertex_shared_global_dof[ii] " << vertex_shared_global_dof[ii]<<std::endl;
       // if(global_to_local_vertex.count(vertex_shared_global_dof[ii])==1)
           {}// do nothing
      //  else
           {}//vertex_color[global_to_local_vertex[vertex_shared_global_dof[ii]]]=vertex_shared_color[ii];
           }
    }

    MPI_Barrier(MPI_GHOST_GROUP);     
}
  
    // MPI_Recv(&vertex_shared_color[0],count_vertex,MPI_INT,world_rank,world_rank,MPI_COMM_WORLD,&status);

    //MPI_Bcast(&mandami,1,MPI_INT,0,MPI_COMM_GHOST);
 
 


MPI_Request request[9];
MPI_Status status[9];
std::vector<int> ricevimi(9,0);
std::vector<int> mandami(9,1);
 
if(world_rank==0)
  all_shared_ranks={0,4,5,8};
if(world_rank==1)
  all_shared_ranks={1,4,6,8};
if(world_rank==2)
  all_shared_ranks={2,6,7,8};
if(world_rank==3)
  all_shared_ranks={3,5,7,8};
if(world_rank==4)
  all_shared_ranks={0,1,4,5,6,8};  
if(world_rank==5)
  all_shared_ranks={0,3,4,5,7,8};
if(world_rank==6)
  all_shared_ranks={1,2,4,6,7,8};
if(world_rank==7)
  all_shared_ranks={2,3,5,6,7,8};
if(world_rank==8)
  all_shared_ranks={0,1,2,3,4,5,6,7,8};
  
std::string s1,s2,s3,output_name;
s1="sending_"; s2=std::to_string(world_rank); s3=".txt";
output_name = s1 + s2 + s3;
std::ofstream outputFileSend(output_name, std::ofstream::out);

if(world_rank<4)
for(int ii=0;ii<4;ii++)
ricevimi[ii]=1;

   for(int ii=0;ii<all_shared_ranks.size();ii++)
     {// receive from processes with a lower rank
      if(all_shared_ranks[ii]<world_rank)
         {MPI_Status status_now;
         int ricevimi_loc;
         int tag=world_rank*world_size+all_shared_ranks[ii];
         MPI_Recv(&ricevimi_loc,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD,&status_now);
         ricevimi[status_now.MPI_SOURCE]=ricevimi_loc;
         //MPI_Irecv(&ricevimi[all_shared_ranks[ii]],1,MPI_INT,all_shared_ranks[ii],tag, MPI_COMM_WORLD,&request[all_shared_ranks[ii]]);
         //MPI_Wait(&request[all_shared_ranks[ii]],NULL);
         //MPI_Recv(&ricevimi_loc,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD,&status[ii]);
         }
      }
    for(int ii=0;ii<all_shared_ranks.size();ii++)  
      if(all_shared_ranks[ii]<world_rank)
         {
          MPI_Wait(&request[all_shared_ranks[ii]],NULL);
          outputFileSend<<"iter=="<<ii<<", world_rank=="<<", "<<world_rank<<" RICEVO "<<(ricevimi[all_shared_ranks[ii]])<<" DA "<<all_shared_ranks[ii]<<"\n";}
         
         
    int sum=0;
      for(int ii=0;ii<all_shared_ranks.size();ii++)
      {
      if(world_rank==5)
       std::cout<<"ricevimi[all_shared_ranks[ii]] "<<ricevimi[all_shared_ranks[ii]]<<std::endl;
      if(all_shared_ranks[ii]<world_rank)  
       {
       sum=sum+ricevimi[all_shared_ranks[ii]];
       
       }
      
      if(world_rank>3 && all_shared_ranks[ii]>world_rank)
      mandami[all_shared_ranks[ii]]=sum;  
      } 
      // coloring
      for(int ii=0;ii<all_shared_ranks.size();ii++)  
      if(all_shared_ranks[ii]==world_rank)
      {
      //std::cout<<"iter=="<<ii<<", world_rank=="<<", "<<world_rank<<", mandami "<<(mandami)<<" ricevimi "<<ricevimi<<std::endl;
      }
      //MPI_Wait(&request, MPI_STATUS_IGNORE);
      // send to processes with a lower rank
      
     // std::cout<<"world_rank=="<<world_rank<<"ricevimi=="<<ricevimi[ii]<<" mandami=="<<mandami<<std::endl;
    for(int ii=0;ii<all_shared_ranks.size();ii++)
      if(all_shared_ranks[ii]>world_rank)
        {
        int tag=world_rank*world_size+all_shared_ranks[ii];
        //MPI_Isend(&mandami[all_shared_ranks[ii]],1,MPI_INT,all_shared_ranks[ii],tag,MPI_COMM_WORLD,&request[all_shared_ranks[ii]]);
        //std::cout<<"SEND Da=="<<world_rank<<" A "<<all_shared_ranks[ii]<<std::endl;
        //MPI_Isend(&mandami[all_shared_ranks[ii]],1,MPI_INT,all_shared_ranks[ii],tag,MPI_COMM_WORLD,&request[all_shared_ranks[ii]]);
        MPI_Send(&mandami[all_shared_ranks[ii]],1,MPI_INT,all_shared_ranks[ii],tag,MPI_COMM_WORLD);
        if(world_rank==4)
         {
         std::cout<<"SEND=="<<mandami[all_shared_ranks[ii]]<<" A=="<<all_shared_ranks[ii]<<" DA "<<world_rank <<std::endl;
         }
        outputFileSend<<"iter=="<<ii<<", world_rank=="<<", "<<world_rank<<", MANDO "<<(mandami[all_shared_ranks[ii]])<<" A "<<all_shared_ranks[ii]<<"\n";
        }


//    {
//     if(world_rank==all_shared_ranks[ii])
//     {
//     for(int jj=ii+1;jj<all_shared_ranks.size();jj++)
//        {
//         MPI_Send(&mandami, 1, MPI_INT, all_shared_ranks[jj], all_shared_ranks[jj], MPI_COMM_WORLD);
//          std::cout<<"iter=="<<ii<<"--IO "<<", "<<world_rank<<", MANDO "<<mandami<<" A "<<jj<<", "<<  prova(jj)<<std::endl;
//        }
//        
//     for(int jj=0;jj<ii;jj++)
//        {
//         MPI_Recv(&ricevimi, 1, MPI_INT,all_shared_ranks[jj], all_shared_ranks[jj], MPI_COMM_WORLD,&status);
//         std::cout<<"--IO "<<", "<<world_rank<<", MANDO "<<mandami<<" A "<<jj<<", "<<  prova(jj)<<std::endl;
//        }       
//        
//        
//     }
// 
//      MPI_Barrier(MPI_COMM_WORLD); 
//     }
    
    
    
    
    
    
       
//   for(int ii=0;ii<world_size-1;ii++)
//     {for(int jj=ii+1;jj<world_size;jj++)
//       if(world_rank==ii)
//        {
//         MPI_Send(&mandami, 1, MPI_INT, jj, jj, MPI_COMM_WORLD);
//          std::cout<<"--IO "<<", "<<world_rank<<", MANDO "<<mandami<<" A "<<jj<<", "<<  prova(jj)<<std::endl;
//        }
//     if(world_rank>ii)
//     {
//     MPI_Recv(&ricevimi, 1, MPI_INT,ii, world_rank, MPI_COMM_WORLD,&status);
//     //std::cout<<all_shared_ranks.size()<<", "<<world_size<<"--IO "<<", "<<world_rank<<", RICEVO "<<ricevimi<<" DA "<<all_shared_ranks[ii]<<std::endl;
//     }
//     
//      MPI_Barrier(MPI_COMM_WORLD); 
//     }
    
    
    
    
    
    
    
    
   //  for(int ii=0;ii<1;ii++)
//     {
//      if(ii!=world_rank) 
//      MPI_Send(&mandami, 1, MPI_INT, world_rank+1, ii, MPI_COMM_WORLD);
//      else
//      MPI_Recv(&ricevimi, 1, MPI_INT, world_rank, ii, MPI_COMM_WORLD,&status);
//          std::cout<<"RICEVO "<<ii<<", "<<world_rank<<", "<<ricevimi<<std::endl;
//      }




s1="coloring_"; s2=std::to_string(world_rank); s3=".txt";
output_name = s1 + s2 + s3;
std::ofstream outputFile(output_name, std::ofstream::out);
for(int ii=0;ii<mesh->num_vertices();ii++)
{
    auto actual_vertex=Vertex(*mesh, ii);
 	auto point_vertex=actual_vertex.point(); 
outputFile << point_vertex[0]<<", "<<point_vertex[1]<<", "<<vertex_color[ii]<<"\n";
    
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

