


#include <dolfin.h>
#include "MixedPoisson.h"
#include "Prova.hpp"
#include <fstream> 
#include "input.hpp"
#include "MeshUniformRefinement.hpp"
#include "Topology.hpp"
#include "RaviartThomasBasisFunctions.hpp"
using namespace dolfin;




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
auto mesh = std::make_shared<UnitSquareMesh>(15,15);
auto gdim = mesh->geometry().dim();  
unsigned int number_of_levels=3;
UniformRefinement uniformrefinement(mesh,number_of_levels,gdim);
auto mesh_list=uniformrefinement.mesh();
auto topology_list=uniformrefinement.topology();
// --The domain is subdivided into world_size subdomains
// 	 each domain has exactly num_domain_vertices vertices
//	 but each processor will contain num_all_vertices
// 	 because also ghost cells are present
// --A vertex belongs to the domain if index_vertex<num_domain_vertices
unsigned int num_domain_vertices=mesh->topology().ghost_offset(0);
unsigned int num_all_vertices=mesh->num_vertices();
auto shared_vertices=mesh->topology().shared_entities(0);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////----------------------------- MPI ----------------------------/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int world_rank,world_size;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  
MPI_Comm_size(MPI_COMM_WORLD, &world_size);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////------------------- TOPOLOGICAL INFORMATION -------------------////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Topology topology(mesh,gdim);
MeshConnectivity topology_N2N= topology.N2N();

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
auto Patch=topology.N2PatchDofs( mesh, dofmap, topology.N2F());

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////--------------------------- COLORING --------------------------////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
unsigned int max_num_colors;    
std::vector< std::vector< unsigned int> > color_2_vertex=color2vertex(mesh,topology.N2N(),shared_vertices,max_num_colors);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////------------------------ SYSTEM ASSEMBLING  -----------------------////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

unsigned int cont=0;
PetscScalar minusone=-1;
PetscScalar one=1;
// this processor owns from first_row to last_row rows of the matrix A
// this processor owns from first_component to last_component rows of the vector b 
PetscReal norm_res;
PetscInt first_row; 
PetscInt last_row; 
PetscInt first_component;
PetscInt last_component;


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

assemble(A,a); //assemble the bilinear form
assemble(b,L); //assemble the external linear form
assemble(bbc,Lzerobc); // bbc=0;
assemble(bbcfound,Lzerobc); // bbcfound=0;

bc.apply(A,b); // apply bcs in a non-symmetric way
bc.apply(bbc); // bbc contains only the bcs
bc1.apply(bbcfound); // bc1 contains now 1 the the bcs dofs position (used to find later the bcs-dofs position)
Mat Amat=A.mat(); 
Vec bvec=b.vec();
Vec bbcvec; 
Vec bbcfoundvec=bbcfound.vec();
Vec Abbcvec;
MatCreateVecs(Amat,&bbcvec,&Abbcvec);



PetscInt globalsize,localsize;
ierr=VecGetSize(bvec,&globalsize);CHKERRQ(ierr);
ierr=VecGetLocalSize(bvec,&localsize);CHKERRQ(ierr);
Vec Ax,res,x; // Ax=A*x, res=b-A*x, x solution, which we initialize to b so that it satisfies bc
VecDuplicate(bvec,&res);
VecDuplicate(bvec,&x);






MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);

bbcvec=bbc.vec();
MatMult(Amat,bbcvec,Abbcvec);
VecAXPY(bvec,minusone,Abbcvec);
VecDestroy(&Abbcvec);







ierr=MatGetOwnershipRange(Amat,&first_row,&last_row);CHKERRQ(ierr);
ierr=VecGetOwnershipRange(res,&first_component,&last_component);CHKERRQ(ierr);


PetscInt L2G_dim=local_to_global_map.size();
Mat* Aloc; // the diagonal submatrix related to all dofs belonging to world_rank + ghost dofs
Vec corrloc,resloc;  // the subvector related to all dofs belonging to world_rank + ghost dofs
VecCreate(PETSC_COMM_SELF,&corrloc);
VecCreateSeq(PETSC_COMM_SELF ,L2G_dim,&corrloc);


IS all_is, patch_is,patch_local_is,all_local_is; 
ierr=createbc(local_to_global_map, Amat,bbcfoundvec,all_is,all_local_is);CHKERRQ(ierr);

Vec vecprova;

ierr=MatCreateSubMatrices(Amat,1,&all_is,&all_is, MAT_INITIAL_MATRIX,&Aloc);CHKERRQ(ierr); // the diagonal submatrix


VecScatter scatter;
ierr=VecScatterCreate(corrloc,NULL,x,all_is,&scatter);CHKERRQ(ierr); 


MatCreateVecs(Amat,&x,&Ax);
VecZeroEntries(x);
MatMult(Amat,x,Ax);
VecWAXPY(res,minusone,Ax,bvec);
ierr=VecGetSubVector(res,all_is,&resloc);CHKERRQ(ierr);


unsigned int smoothing_steps=2000;

ierr=ColoredGaussSeidel(Amat,x,bvec,Ax,res,smoothing_steps,max_num_colors,color_2_vertex,Patch,corrloc,resloc,Aloc,all_local_is,all_is);CHKERRQ(ierr); 

ierr=VecScatterDestroy(&scatter);CHKERRQ(ierr); 





VecAXPY(x,one,bbcvec);


//Mat WC2WF=WC2WFMatrix2D(mesh_list[1],mesh_list[2],topology_list[1],topology_list[2],Amat);



Vec xsol;
PetscScalar *arraylocalghosted;
ierr=VectorToGhostedVector(x,xsol,arraylocalghosted,all_is,local_to_global_map,L2G_dim,globalsize,localsize );CHKERRQ(ierr);
auto sol = std::make_shared<PETScVector>(xsol); 
Function w(W,sol);
// Extract sub functions (function views)
Function& sigma = w[0];
Function& u = w[1];
 
File fileu("../output/mixedpoisson_u.pvd");
fileu << u;
File filesigma("../output/mixedpoisson_sigma.pvd");
filesigma << sigma;
 
//  auto W2=std::make_shared<MixedPoisson::FunctionSpace>(mesh_list[2]);
// Function w2(W2);
// Function& sigma2 = w2[0];
//    File file2("../output/mixedpoisson_prova2.pvd");
//    file2 << sigma2;
//    
//  auto W1=std::make_shared<MixedPoisson::FunctionSpace>(mesh_list[1]);
// Function w1(W1);
// Function& sigma1 = w1[0];
//    File file1("../output/mixedpoisson_prova1.pvd");
//    file1 << sigma1;
// 
// Function w0(W);
// Function& sigma0 = w0[0];
//    File file0("../output/mixedpoisson_prova0.pvd");
//    file0 << sigma0;      


PetscFree(arraylocalghosted);
VecDestroy(&xsol);
ISDestroy(&all_is);
PetscFree(all_is);
VecDestroy(&Ax);
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





//ierr=ISView(all_is,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
//ierr=VecView(resloc, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector


//WriteMatrixToFile("Aasym.txt","Aasym",W_dim,W_dim,Amat,(world_rank==0));
//WriteVectorToFile("basym.txt","basym",W_dim,bvec,(world_rank==0));

//ierr=VecView(bvec, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector
//ierr=VecView(Abbcvec, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector
//ierr=VecView(bbcfoundvec, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr); // view the subvector
//WriteVectorToFile("indices"+std::to_string(world_rank)+".txt","indices"+std::to_string(world_rank),L2G_dim,vecprova,1);





// Function w(W);
// solve(a == L, w, bc);