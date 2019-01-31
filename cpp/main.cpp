


#include <dolfin.h>
#include "MixedPoisson.h"
#include "Prova.hpp"
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

class ZeroSource : public Expression
{
  public: 
  ZeroSource(): Expression(2) {}
  void eval(Array<double>& values, const Array<double>& x) const
  {

    values[0] = 0.0;
    values[1] = 0.0;
  }
};




// Boundary source for flux boundary condition
class BoundarySource : public Expression
{
public:

  BoundarySource(const Mesh& mesh) : Expression(2), mesh(mesh) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 1.0;
    values[1] = 1.0;
  }

private:

  const Mesh& mesh;

};

class BoundarySource1 : public Expression
{
public:

  BoundarySource1(const Mesh& mesh) : Expression(2), mesh(mesh) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 1.0;
    values[1] = 1.0;
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
//////////////////////////////////////////
/////////     CONFIGURATION       ////////
//////////////////////////////////////////
parameters["ghost_mode"] = "shared_vertex";
parameters["reorder_vertices_gps"] = true;
parameters["reorder_cells_gps"] = true;
parameters["linear_algebra_backend"] = "PETSc";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////----------------------------- MESH ---------------------------/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
auto mesh = std::make_shared<UnitSquareMesh>(5,5);
// --The domain is subdivided into world_size subdomains
// 	 each domain has exactly num_domain_vertices vertices
//	 but each processor will contain num_all_vertices
// 	 because also ghost cells are present
// --A vertex belongs to the domain if index_vertex<num_domain_vertices
unsigned int num_domain_vertices=mesh->topology().ghost_offset(0);
unsigned int num_all_vertices=mesh->num_vertices();
auto shared_vertices=mesh->topology().shared_entities(0);
auto gdim = mesh->geometry().dim();  

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////----------------------------- MPI ----------------------------/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int world_rank,world_size;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
MPI_Comm_size(MPI_COMM_WORLD, &world_size);









    
    
    std::vector<unsigned int> shared_rank, all_shared_ranks;
    std::map<unsigned int, unsigned int> shared_rank_map;
    
    find_communicating_processes(shared_rank,all_shared_ranks, shared_vertices,shared_rank_map);
    std::map<unsigned int, unsigned int> global_to_local_vertex=global2local_vertex(mesh,shared_vertices) ;

  
  std::vector<bool> use_vertex_color=divide_nodes_among_processes(mesh,shared_vertices,shared_rank,all_shared_ranks, shared_rank_map,global_to_local_vertex);


   
    for(int jj=0;jj<world_size;jj++)
    {
    if(world_rank==jj)
    {
        for(int mm=0;mm<use_vertex_color.size();mm++)
        std::cout<<" rank=="<<world_rank<<", use_vertex_color["<<mm<<"]="<<use_vertex_color[mm]<<std::endl;
    }
        std::cout<<std::endl;
       MPI_Barrier(MPI_COMM_WORLD);
       }
    


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////------------------- TOPOLOGICAL INFORMATION -------------------////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // node to edge 
mesh->init(0,gdim-1);
MeshConnectivity topology_N2F=mesh->topology()(0,gdim-1);
// face to node
mesh->init(gdim-1,0);
MeshConnectivity topology_F2N=mesh->topology()(gdim-1,0);
// element to node
mesh->init(gdim,0);
MeshConnectivity topology_K2N=mesh->topology()(gdim,0);
// element to face
mesh->init(gdim,gdim-1);
MeshConnectivity topology_K2F=mesh->topology()(gdim,gdim-1);
// face to element
mesh->init(gdim-1,gdim);
MeshConnectivity topology_F2K=mesh->topology()(gdim-1,gdim);
// node to element
mesh->init(0,gdim);
MeshConnectivity topology_N2K=mesh->topology()(0,gdim);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////------------------------ FUNCTION SPACE  ----------------------////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

auto W = std::make_shared<MixedPoisson::FunctionSpace>(mesh);
PetscInt W_dim=W->dim(); 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////---------------------------- DOFMAP ---------------------------////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<const dolfin::GenericDofMap> dofmap = W->dofmap();
std::vector<std::size_t> local_to_global_map;
dofmap->tabulate_local_to_global_dofs(local_to_global_map);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////---------------------- TOPOLOGICAL AND DOFMAP INFORMATION -------------------//////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// map between a local vertex and the local vertices belonging which share an ege with it
std::vector<std::vector< int > > topology_N2N=topologyN2N(mesh, topology_N2F, topology_F2N);
// map between a local vertex and the local dofs belonging to its patch
std::vector<std::vector< unsigned int > > Patch=topologyN2PatchDofs( mesh, dofmap, topology_N2F);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////--------------------------- COLORING --------------------------////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
unsigned int max_num_colors;    
std::vector< std::vector< unsigned int> > color_2_vertex=color2vertex(mesh,topology_N2N,shared_vertices,max_num_colors);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////------------------------ OUTPUT INFORMATION -----------------------////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::ofstream myfile;
std::string  Amatfiletxt="Amatfile.txt",resfile="resfile",resfiletxt=resfile+".txt";
std::string Alocfile="Alocfile",resfileloc="resfileloc", Apatchfile="Apatchfile";
std::string txt=".txt";
std::string stringrank=std::to_string(world_rank);
Alocfile=Alocfile+stringrank;
Apatchfile=Apatchfile+stringrank+"_";
resfileloc=resfileloc+stringrank;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////------------------------ SYSTEM ASSEMBLING  -----------------------////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int cont=0;
PetscScalar minusone=-1;
const double Lambda = 1.0;
const double Mu = 1.0;
const double Beta = 1.0/(2*Mu);

auto C_equilibrium=std::make_shared<Constant>(1.0);
auto C_zero=std::make_shared<Constant>(0.0);
auto C_constitutive=std::make_shared<Constant>(1.0);
auto beta=std::make_shared<Constant>(1.0/(2*Mu));
auto alpha=std::make_shared<Constant>(((-Beta) * Lambda /(2*Lambda +2*Mu)));
 
MixedPoisson::BilinearForm a(W, W);
MixedPoisson::LinearForm L(W);
MixedPoisson::LinearForm Lzerobc(W);



a.C_constitutive=C_constitutive;
a.C_equilibrium=C_equilibrium;
a.beta=beta;
a.alpha=alpha;

auto f = std::make_shared<Source>();
auto fzero = std::make_shared<ZeroSource>();

L.f = f;
L.C_equilibrium=C_equilibrium;
Lzerobc.f=fzero;
Lzerobc.C_equilibrium=C_zero;
// Define boundary condition
auto G = std::make_shared<BoundarySource>(*mesh);
auto G1 = std::make_shared<BoundarySource1>(*mesh);
auto boundary = std::make_shared<EssentialBoundary>();
DirichletBC bc(W->sub(1), G, boundary);
DirichletBC bc1(W->sub(1), G1, boundary);


 








PetscErrorCode ierr; 
PETScMatrix A; // define parallel matrix A
PETScVector b,bbc,bbcfound; // define parallel vector b

assemble(A,a);
assemble(b,L);
assemble(bbc,Lzerobc);
assemble(bbcfound,Lzerobc);

bc.apply(A,b);
bc.apply(bbc);
bc1.apply(bbcfound);
Mat Amat=A.mat(); 
Vec bvec=b.vec();
Vec bbcvec;
Vec bbcfoundvec=bbcfound.vec();
Vec Abbcvec;
MatCreateVecs(Amat,&bbcvec,&Abbcvec);

WriteMatrixToFile("Aasym.txt","Aasym",W_dim,W_dim,Amat,(world_rank==0));
WriteVectorToFile("basym.txt","basym",W_dim,bvec,(world_rank==0));


PetscInt globalsize;
VecGetSize(bvec,&globalsize);
Vec Ax,res,x; // Ax=A*x, res=b-A*x, x solution, which we initialize to b so that it satisfies bc

VecDuplicate(bvec,&res);
VecDuplicate(bvec,&x);



// this processor owns from first_row to last_row rows of the matrix A
// this processor owns from first_component to last_component rows of the vector b 
PetscReal norm_res;
PetscInt first_row; 
PetscInt last_row; 
PetscInt first_component;
PetscInt last_component;


MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
//ierr=MatView(Amat, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
//ierr=VecView(bvec, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);

bbcvec=bbc.vec();
MatMult(Amat,bbcvec,Abbcvec);
VecAXPY(bvec,minusone,Abbcvec);




//ierr=VecView(bvec, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector
//ierr=VecView(Abbcvec, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector
//ierr=VecView(bbcfoundvec, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector




ierr=MatGetOwnershipRange(Amat,&first_row,&last_row);CHKERRQ(ierr);
ierr=VecGetOwnershipRange(res,&first_component,&last_component);CHKERRQ(ierr);




PetscInt localsize=last_row-first_row;//-1; <-------------------- check here
PetscInt L2G_dim=local_to_global_map.size();

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
//std::cout<<"world_rank=="<<world_rank<<", localsize=="<<localsize<<", bc_local_dofs="<<bc_local_dofs[ii]<<"bc_dofs="<<std::endl;
}
PetscInt* bc_dofs=(PetscInt*)malloc(sizeof(PetscInt)*cont);

cont=0;
for(int ii=0;ii<localsize;ii++)
{if(bc_local_dofs[ii]==1)
    {bc_dofs[cont]=one_to_localsize_vec[ii];
     cont++;} 
}

for(int ii=0;ii<cont;ii++)
{}//std::cout<<"world_rank=="<<world_rank<<", bc_dofs="<<bc_dofs[ii]<<" bc_dofs_size="<<cont<<std::endl;
MatZeroRows(Amat,cont,bc_dofs,1,0,0);
MatZeroRowsColumns(Amat,cont,bc_dofs,1,0,0);
//ierr=MatView(Amat, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); 
cont=0;

Mat* Aloc; // the diagonal submatrix related to all dofs belonging to world_rank + ghost dofs
Vec corrloc,resloc;  // the subvector related to all dofs belonging to world_rank + ghost dofs
VecCreate(PETSC_COMM_SELF,&corrloc);
VecCreateSeq(PETSC_COMM_SELF ,L2G_dim,&corrloc);


IS all_is, patch_is,patch_local_is,all_local_is; // all_is=IS related to all dofs belonging to world_rank + ghost dofs and patch_is=IS related to the patch
ierr=ISCreateGeneral(PETSC_COMM_SELF,L2G_dim,all_indices,PETSC_COPY_VALUES,&all_is);CHKERRQ(ierr);
ierr=ISCreateGeneral(PETSC_COMM_SELF,L2G_dim,one_to_L2G_dim_vec,PETSC_COPY_VALUES,&all_local_is);CHKERRQ(ierr);


Vec vecprova;
VecCreateSeq(PETSC_COMM_SELF,L2G_dim,&vecprova);
PetscScalar* all_indices_scalar=(PetscScalar*)malloc(sizeof(PetscScalar)*L2G_dim);
for(int ii=0;ii<L2G_dim;ii++)
   { 
    all_indices_scalar[ii]=(PetscScalar)all_indices[ii];
    //std::cout<<"all_indices"<<all_indices_scalar[ii]<<std::endl;
    }
    
VecSetValues(vecprova,L2G_dim,one_to_L2G_dim_vec,all_indices_scalar,INSERT_VALUES);    
WriteVectorToFile("indices"+std::to_string(world_rank)+".txt","indices"+std::to_string(world_rank),L2G_dim,vecprova,1);
VecDestroy(&vecprova);



WritetopologyN2PatchDofs(mesh, dofmap, topology_N2F,world_rank,shared_vertices,all_indices_scalar);

//ierr=ISView(all_is,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
//ierr=VecView(resloc, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector

ierr=MatCreateSubMatrices(Amat,1,&all_is,&all_is, MAT_INITIAL_MATRIX,&Aloc);CHKERRQ(ierr); // the diagonal submatrix
VecScatter scatter;

ierr=VecScatterCreate(corrloc,NULL,x,all_is,&scatter);CHKERRQ(ierr); 


MatCreateVecs(Amat,&x,&Ax);
VecZeroEntries(x);
MatMult(Amat,x,Ax);
VecWAXPY(res,minusone,Ax,bvec);

std::string matrixname="Amat";
WriteMatrixToFile(Amatfiletxt,matrixname,W_dim,W_dim,Amat,(world_rank==0));
WriteVectorToFile(resfiletxt,resfile,W_dim,bvec,(world_rank==0));
WritePatchToFile("Patch",Patch,world_rank,max_num_colors,color_2_vertex,all_indices_scalar);
VecNorm(res, NORM_2,&norm_res);


std::cout<<"world_rank="<<world_rank<<", num_domain_vertices="<<num_domain_vertices<<std::endl;
std::cout<<"world_rank="<<world_rank<<", num_all_vertices="<<num_all_vertices<<std::endl;
std::cout<<"world_rank="<<world_rank<<", mesh->topology().ghost_offset(0)="<<mesh->topology().ghost_offset(0)<<std::endl;




if(world_rank==0)
{std::cout.precision(14);
std::cout<<std::scientific<<"smoothing_step=="<<0<<", norm_res=="<<norm_res<<std::endl;}

//WriteMatrixToFile(Alocfile+txt,Alocfile,L2G_dim,L2G_dim,*Aloc,1);

ierr=VecGetSubVector(res,all_is,&resloc);CHKERRQ(ierr);
//WriteVectorToFile(resfileloc+txt,resfileloc,L2G_dim,resloc,1);



PetscFree(all_indices_scalar);
PetscFree(one_to_L2G_dim_vec);
PetscFree(all_indices);
unsigned int smoothing_steps=500;
for(int ss=0;ss<smoothing_steps;ss++)
for(int cc=1;cc<max_num_colors;cc++)
{
// extract the local residual and set to zero the local correction
ierr=VecGetSubVector(res,all_is,&resloc);CHKERRQ(ierr); // the subvector
VecZeroEntries(corrloc);
//loop on internal nodes

for(std::vector<unsigned int>::iterator iter_c=color_2_vertex[0].begin();iter_c!=color_2_vertex[0].end();iter_c++)
   {
   
   int ii=*iter_c;
   if(world_rank==1)
   std::cout<<"ss="<<ss<<" rank="<<world_rank<< " , color="<<0<<" , "<<ii<<" ,Patch[ii].size()="<<Patch[ii].size()<<std::endl;
        unsigned int vertex_index=ii;
        auto actual_vertex=Vertex(*mesh, vertex_index);
        auto point=actual_vertex.point();
        unsigned int global_vertex_index=actual_vertex.global_index();
   //std::cout<<"rank=="<<world_rank<<", loc="<<ii<<", glob=="<<global_vertex_index<<std::endl;
   
   
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


//loop on boundary nodes (colored)
for(std::vector<unsigned int>::iterator iter_c=color_2_vertex[cc].begin();iter_c!=color_2_vertex[cc].end();iter_c++)
   {
 
   int ii=*iter_c;
   if(world_rank==1)
   std::cout<<"ss="<<ss<<" rank="<<world_rank<< " , color="<<cc<<" , "<<ii<<" ,Patch[ii].size()="<<Patch[ii].size()<<std::endl;
   //std::cout<<" ii= "<<ii<<std::endl;
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
   
   //WriteVectorToFile("corrloc1prima","corrloc1prima",L2G_dim,corrloc,(world_rank==0));
   MatMult(Apatchall,corrloc,Acorrloc);
   //WriteMatrixToFile("Apatchall1","Apatchall1",sizepatch,L2G_dim,Apatchall,1);
  // WriteVectorToFile("respatch1prima","respatch1prima",sizepatch,respatch,(world_rank==0));
   VecAXPY(respatch,minusone,Acorrloc);
   
  // WriteMatrixToFile(Apatchfile+std::to_string(ii)+txt,Alocfile+std::to_string(ii),sizepatch,sizepatch,Apatch,1);
  // WriteVectorToFile("Acorrloc1","Acorrloc1",sizepatch,Acorrloc,(world_rank==0));
  // WriteVectorToFile("respatch1","respatch1",sizepatch,respatch,(world_rank==0));

   //respatch=
   
   //ierr=MatView(Apatch, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the patch submatrix
   //ierr=VecView(respatch, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the patch residual
   //VecGetSize(respatch,&sizepatch);
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
MPI_Barrier(MPI_COMM_WORLD);
 // view the subvector
//VecView(corrloc, PETSC_VIEWER_STDOUT_SELF);
VecNorm(corrloc, NORM_2,&norm_res);
//std::cout<<"smoothing_step=="<<ss<<", color="<<cc<<", norm_corrloc=="<<norm_res<<std::endl;
VecScatterBegin(scatter,corrloc,x, ADD_VALUES,SCATTER_FORWARD);
VecScatterEnd(scatter,corrloc,x, ADD_VALUES,SCATTER_FORWARD);
MatMult(Amat,x,Ax);
VecWAXPY(res,minusone,Ax,bvec);
PetscReal norm_corrloc;
VecNorm(corrloc, NORM_2,&norm_corrloc);
VecNorm(res, NORM_2,&norm_res);
//WriteVectorToFile("corrloc"+std::to_string(world_rank)+"_"+std::to_string(cc),"corrloc"+std::to_string(world_rank)+"_"+std::to_string(cc),L2G_dim,corrloc,1);
//WriteVectorToFile("x_"+std::to_string(cc),"x_"+std::to_string(cc),W_dim,x,(world_rank==0));

if(world_rank==0)
{std::cout.precision(14);
std::cout<<std::scientific<<"smoothing_step=="<<ss<<", color="<<cc<<" norm_corrloc=="<<norm_corrloc<<", norm_res=="<<norm_res<<std::endl;}
VecNorm(x, NORM_2,&norm_res);
if(world_rank==0)
{}//std::cout<<"smoothing_step=="<<ss<<", color="<<cc<<", norm_x=="<<norm_res<<std::endl;
}


//savemat("x.mat", x);//dict(x=x->array(x)));
//PetscViewer  PETSC_VIEWER_MATLAB_(	PETSC_VIEWER_ASCII_MATLAB);
//ierr=VecView(x, PETSC_VIEWER_MATLAB_(PETSC_COMM_WORLD));CHKERRQ(ierr);
ierr=VecScatterDestroy(&scatter);CHKERRQ(ierr); 







  
	//std::cout<<"rank=="<<world_rank<<" ss="<<jj<<" cc=="<<cc<<"/"<<used_color.size()<<std::endl;
	//MPI_Barrier(MPI_COMM_WORLD);
	//PetscScalar one=1;
	 //VecAXPY(x,one,Vec x);
	//}

ISDestroy(&all_is);
PetscFree(all_is);





















  // Compute solution
Function w(W);//,x);
 // copyx);
 // solve(a == L, w, bc);

// Now, the separate components ``sigma`` and ``u`` of the solution can
// be extracted by taking components. These can then be written to XDMF
// files.
// 
// .. code-reslock:: cpp

  // Extract sub functions (function views)
   Function& sigma = w[0];
   Function& u = w[1];
 
   //File file1("mixedpoisson_u.pvd");
   //file1 << u;
   //File file2("mixedpoisson_sigma.pvd");
   //file2 << sigma;
  
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




/* 
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
	 

 */

//   for(int jj=0;jj<world_size;jj++)
//      {
//      MPI_Barrier(MPI_COMM_WORLD);
//      if(world_rank==jj)
// while(it != shared_cells.end())
// {
// auto nodes=topology_K2N(it->first);
// //std::cout<<"world_rank: "<<world_rank<<", cell "<<it->first<<", nodes "<< nodes[0]<<", "<<nodes[1]<<", "<<nodes[2]<<std::endl;
// for(int ii=0;ii<3;ii++)
// {
// auto actual_vertex=Vertex(*mesh, nodes[ii]);
// auto point_vertex=actual_vertex.point();
// //std::cout<<point_vertex[0]<<", "<<point_vertex[1]<<std::endl;
// }
// it++;
// }
// }


// MatView(submat, PETSC_VIEWER_STDOUT_WORLD);

//assemble(b,L);
































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















// Vec vecghost;
// VecCreateGhost(PETSC_COMM_WORLD,localsize,globalsize,nghost,ghost_indices,&vecghost);
// ierr=VecSetValues(vecghost,nghost,ghost_indices,ghost_scalars,ADD_VALUES);CHKERRQ(ierr);
// VecGhostUpdateBegin(vecghost,ADD_VALUES,SCATTER_REVERSE);
// VecGhostUpdateEnd(vecghost,ADD_VALUES,SCATTER_REVERSE);
// VecGhostUpdateBegin(vecghost,INSERT_VALUES,SCATTER_FORWARD);
// VecGhostUpdateEnd(vecghost,INSERT_VALUES,SCATTER_FORWARD);
// VecAssemblyBegin(vecghost);
// VecAssemblyEnd(vecghost);


// VecGhostUpdateBegin(vecghost, ADD_VALUES, SCATTER_FORWARD);
// VecGhostUpdateEnd(vecghost, ADD_VALUES, SCATTER_FORWARD);
//
// 
// //VecCopy(bvec,vecghost);
//ierr=VecView(vecghost, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector



//PetscScalar *v;
//PetscScalar *z;
//PetscInt *indices;
//PetscInt n_loc=4;
//indices=(PetscInt *)malloc(n_loc*n_loc*sizeof(PetscInt)); 
//PetscMalloc1(n_loc,&indices); //allocate the space for vector 
//PetscMalloc1(n_loc*n_loc,&v); //allocate the space for matrix 
//PetscMalloc1(n_loc,&z); //allocate the space for vector 























// std::unordered_map<long unsigned int, double> boundary_dofs_and_values;
// bc.get_boundary_values(boundary_dofs_and_values);
// unsigned int boundary_dofs_and_values_size=boundary_dofs_and_values.size();
// PetscInt* boundary_dofs=(PetscInt* )malloc(sizeof(PetscInt)*boundary_dofs_and_values_size);
// PetscScalar* boundary_values=(PetscScalar*)malloc(sizeof(PetscScalar)*boundary_dofs_and_values_size);
// 
// 
// std::cout << "rank="<<world_rank<<"boundary_values size="<<boundary_dofs_and_values_size << std::endl;
// for(std::pair<long unsigned int, double> element : boundary_dofs_and_values)
// {
// boundary_dofs[cont]=element.first;
// boundary_values[cont]=element.second;
// cont++;
// }
// 
// for(int ii=0;ii<cont;ii++)
// std::cout << "rank="<<world_rank<<"boundary_dofs="<<boundary_dofs[ii] << ", boundary_values=" << boundary_values[ii] << std::endl;
// 
// cont=0;
// 
// std::ofstream doffile;
//  doffile.open("dof"+std::to_string(world_rank)+".txt");
//  for(int cc=0;cc<max_num_colors;cc++)  
//     {doffile <<"dof{"<<std::to_string(world_rank+1)+"}{";doffile<<std::to_string(cc+1);doffile<<"}=[";
//    for(int ii=0;ii<color_2_vertex[cc].size();ii++)  
//         {
//         unsigned int vertex_index=color_2_vertex[cc][ii];
//         auto actual_vertex=Vertex(*mesh, vertex_index);
//         auto point=actual_vertex.point();
//         unsigned int global_vertex_index=actual_vertex.global_index();
//         doffile <<std::to_string(point[0]);doffile <<", "; doffile <<std::to_string(point[1]);doffile <<", ";doffile << cc; doffile <<", ";doffile << std::to_string(global_vertex_index);doffile <<";\n";
//         }
//         doffile<<"];\n";
//      }
// 
// doffile.close();
// PetscFree(boundary_dofs);
// PetscFree(boundary_values);









//auto ke=boundary_values->keys();

// Create(PETSC_COMM_WORLD,&res);
// VecSetSizes(res,PETSC_DECIDE, globalsize);
// VecSetUp(res);

// VecCreate(PETSC_COMM_WORLD,&x);
// VecSetSizes(x,PETSC_DECIDE, globalsize);
// VecSetUp(x);




//std::pair<long unsigned int, long unsigned int> ownership_range=dofmap->ownership_range();

 //for(int ii=0;ii<local_to_global_map.size();ii++)
 //{
 //std::cout<<"---xxxx------xxxx------xxxx------xxxx---- world_rank "<<world_rank<<", component="<<ii<<", local_to_global_map="<<local_to_global_map[ii]<<std::endl;
 //}

//std::cout<<"--------------->world_rank "<<world_rank<<"W->dim()=="<<W->dim()<<", mesh->num_vertices() "<<mesh->num_vertices()<<" dofmap->ownership_range=["<<ownership_range.first<<", "<<ownership_range.second<<"] "<<std::endl;

//std::map<int, std::set<unsigned int>>::iterator it_sharedvertex = shared_vertices.begin();
