// Mixed formulation for Poisson equation (C++)
// ============================================
// 
// This demo illustrates how to solve Poisson equation using a mixed
// (two-field) formulation. In particular, it illustrates how to
// 
// * Use mixed and non-continuous finite element spaces
// * Set essential boundary conditions for subspaces and H(div) spaces
// * Define a (vector-valued) expression using additional geometry information
// 
// 
// Equation and problem definition
// -------------------------------
// 
// An alternative formulation of Poisson equation can be formulated by
// introducing an additional (vector) variable, namely the (negative)
// flux: :math:`\sigma = \nabla u`. The partial differential equations
// then read
// 
// .. math::
//    \sigma - \nabla u &= 0 \quad {\rm in} \ \Omega, \\
//    \nabla \cdot \sigma &= - f \quad {\rm in} \ \Omega,
// 
// with boundary conditions
// 
// .. math::
//    u = u_0 \quad {\rm on} \ \Gamma_{D},  \\
//    \sigma \cdot n = g \quad {\rm on} \ \Gamma_{N}.
// 
// The same equations arise in connection with flow in porous media, and
// are also referred to as Darcy flow.
// 
// After multiplying by test functions :math:`\tau` and :math:`v`,
// integrating over the domain, and integrating the gradient term by
// parts, one obtains the following variational formulation: find
// :math:`\sigma \in \Sigma` and :math:`v \in V` satisfying
// 
// .. math::
//    \int_{\Omega} (\sigma \cdot \tau + \nabla \cdot \tau \ u) \ {\rm d} x
//    &= \int_{\Gamma} \tau \cdot n \ u \ {\rm d} s
//    \quad \forall \ \tau \in \Sigma, \\
// 
//    \int_{\Omega} \nabla \cdot \sigma v \ {\rm d} x
//    &= - \int_{\Omega} f \ v \ {\rm d} x
//    \quad \forall \ v \in V.
// 
// Here :math:`n` denotes the outward pointing normal vector on the
// boundary. Looking at the variational form, we see that the boundary
// condition for the flux (:math:`\sigma \cdot n = g`) is now an
// essential boundary condition (which should be enforced in the function
// space), while the other boundary condition (:math:`u = u_0`) is a
// natural boundary condition (which should be applied to the variational
// form). Inserting the boundary conditions, this variational problem can
// be phrased in the general form: find :math:`(\sigma, u) \in \Sigma_g
// \times V` such that
// 
// .. math::
// 
//    a((\sigma, u), (\tau, v)) = L((\tau, v))
//    \quad \forall \ (\tau, v) \in \Sigma_0 \times V
// 
// where the variational forms :math:`a` and :math:`L` are defined as
// 
// .. math::
// 
//    a((\sigma, u), (\tau, v)) &=
//      \int_{\Omega} \sigma \cdot \tau + \nabla \cdot \tau \ u
//    + \nabla \cdot \sigma \ v \ {\rm d} x \\
//    L((\tau, v)) &= - \int_{\Omega} f v \ {\rm d} x
//    + \int_{\Gamma_D} u_0 \tau \cdot n  \ {\rm d} s
// 
// and :math:`\Sigma_g = \{ \tau \in H({\rm div}) \text{ such that } \tau \cdot n|_{\Gamma_N} = g \}`
// and :math:`V = L^2(\Omega)`.
// 
// To discretize the above formulation, two discrete function spaces
// :math:`\Sigma_h \subset \Sigma` and :math:`V_h \subset V` are needed
// to form a mixed function space :math:`\Sigma_h \times V_h`. A stable
// choice of finite element spaces is to let :math:`\Sigma_h` be the
// Brezzi-Douglas-Marini elements of polynomial order :math:`k` and let
// :math:`V_h` be discontinuous elements of polynomial order :math:`k-1`.
// 
// We will use the same definitions of functions and boundaries as in the
// demo for Poisson's equation. These are:
// 
// * :math:`\Omega = [0,1] \times [0,1]` (a unit square)
// * :math:`\Gamma_{D} = \{(0, y) \cup (1, y) \in \partial \Omega\}`
// * :math:`\Gamma_{N} = \{(x, 0) \cup (x, 1) \in \partial \Omega\}`
// * :math:`u_0 = 0`
// * :math:`g = \sin(5x)`   (flux)
// * :math:`f = 10\exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)`   (source term)
// 
// With the above input the solution for :math:`u` and :math:`\sigma` will look as
// follows:
// 
// .. image:: ../mixed-poisson_u.png
//     :scale: 75
//     :align: center
// 
// .. image:: ../mixed-poisson_sigma.png
//     :scale: 75
//     :align: center
// 
// 
// Implementation
// --------------
// 
// The implementation is split in two files, a form file containing the definition
// of the variational forms expressed in UFL and the solver which is implemented
// in a C++ file.
// 
// Running this demo requires the files: :download:`main.cpp`,
// :download:`MixedPoisson.ufl` and :download:`CMakeLists.txt`.
// 
// 
// UFL form file
// ^^^^^^^^^^^^^
// 
// The UFL file is implemented in :download:`MixedPoisson.ufl`, and the
// explanation of the UFL file can be found at :doc:`here <MixedPoisson.ufl>`.
// 
// 
// C++ program
// ^^^^^^^^^^^
// 
// The solver is implemented in the :download:`main.cpp` file.
// 
// At the top we include the DOLFIN header file and the generated header
// file containing the variational forms.  For convenience we also
// include the DOLFIN namespace.
// 
// .. code-block:: cpp


#include <dolfin.h>
#include "MixedPoisson.h"


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

auto mesh = std::make_shared<UnitSquareMesh>(2,2);

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
auto colors=mesh->color("vertex");

if(world_rank==1)
{
for(int ii=0;ii<colors.size();ii++)
{auto cell=Cell(*mesh,ii);
std::vector<double> cell_coor;
cell.get_vertex_coordinates(cell_coor);
for(int jj=0;jj<0.5*cell_coor.size();jj++)
{
//std::cout<<"world_rank["<<world_rank<<"]: "<<"colors["<<ii<<"]: "<<colors[ii]<<" cell "<< cell_coor[2*jj]<<", "<< cell_coor[2*jj+1]<< std::endl;
}

}


}
 
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



std::vector<unsigned int> keys;
std::vector<std::set<unsigned int> > vals;
    for (std::map< int, std::set<unsigned int> >::iterator it=shared_vertices.begin(); it!= shared_vertices.end(); ++it)
       {keys.push_back(it->first);
       vals.push_back(it->second);
       }




   
  for(int jj=0;jj<world_size;jj++)
     {
      MPI_Barrier(MPI_COMM_WORLD);
     if(world_rank==jj)
          {
          for(int ii=0;ii<keys.size();ii++)
			{
 			 std::cout<<"KEYS NUMBER: "<<ii<<std::endl;
 			 auto vertex_index=keys[ii];
 			 auto actual_vertex=Vertex(*mesh, vertex_index);
 			 auto point_vertex=actual_vertex.point();

 	       std::cout<<"world_rank: "<<world_rank<<"keys "<<keys[ii]<<"  point=(" <<point_vertex[0]<<", "<<point_vertex[1]<<") "<<std::endl;
        	for ( std::set<unsigned int>::iterator setit = vals[ii].begin(); setit != vals[ii].end(); ++setit)
             std::cout<<"--- "<<*setit<<std::endl; 
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
std::cout<<"world_rank: "<<world_rank<<", cell "<<it->first<<", nodes "<< nodes[0]<<", "<<nodes[1]<<", "<<nodes[2]<<std::endl;
for(int ii=0;ii<3;ii++)
{
auto actual_vertex=Vertex(*mesh, nodes[ii]);
auto point_vertex=actual_vertex.point();
std::cout<<point_vertex[0]<<", "<<point_vertex[1]<<std::endl;
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

