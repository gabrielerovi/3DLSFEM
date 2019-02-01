


#include <dolfin.h>
#include <fstream> 
using namespace dolfin;



void WritePatchToFile(const std::string &filename,const std::vector<std::vector< unsigned int > > &Patch, const unsigned int rank, 
                      const unsigned int &max_num_colors,const std::vector< std::vector< unsigned int> > &color_2_vertex,
                      PetscScalar * all_indices_scalar)
{
std::ofstream myfile;
std::string outputname="../output/"+filename+std::to_string(rank);
myfile.open (outputname);
for(int ii=0;ii<Patch.size();ii++)
   {myfile <<"Patch{"+std::to_string(rank+1)+"}{";myfile<<std::to_string(ii+1);myfile<<"}=[";
    for(int jj=0;jj<Patch[ii].size();jj ++)
       {myfile << std::to_string(all_indices_scalar[Patch[ii][jj]]+1);myfile <<" ";}
    myfile <<"];\n";
       }
myfile <<" \n\n";
myfile.close();

myfile.open ("../output/color2vertex"+std::to_string(rank)+".txt");
for(int cc=0; cc<max_num_colors;cc++)
    {myfile <<"color2vertex{"+std::to_string(rank+1)+"}{";myfile<<std::to_string(cc+1);myfile<<"}=[";
     for(int ii=0; ii<color_2_vertex[cc].size();ii++)
        {int vertex=color_2_vertex[cc][ii];
        myfile << std::to_string(vertex+1);myfile <<" ";
        }
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
myfile.open ("../output/"+filename);
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
{N_indices[ii]=ii;}
Vec vecloc;
IS is; // all_is=IS related to all dofs belonging to world_rank + ghost dofs and patch_is=IS related to the patch
ISCreateGeneral(PETSC_COMM_SELF,N,N_indices,PETSC_COPY_VALUES,&is);//CHKERRQ(ierr);
VecGetSubVector(vec,is,&vecloc);
VecGetValues(vecloc,N,N_indices,vectortowrite);
if(right_rank)
{
myfile.open ("../output/"+filename);
myfile <<vectorname; 
myfile <<"=[";
for(int ii=0;ii<N;ii++)
       {myfile << std::to_string(vectortowrite[ii]);myfile <<";\n";}
    myfile <<"\n";
myfile <<"]";
myfile.close();
}
}


void WriteColorToFile(std::shared_ptr<dolfin::Mesh> mesh, const std::vector<unsigned int> &vertex_color,const std::vector<bool> &use_vertex_color)
{
int world_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

int cont=0;
VertexIterator vertexiterator=VertexIterator(*mesh);
for (; !vertexiterator.end(); ++vertexiterator)
 cont++;


std::string s1,s2,s3,output_name;
s1="../output/coloring_"; s2=std::to_string(world_rank); s3=".txt";
output_name = s1 + s2 + s3;
std::ofstream outputFile(output_name, std::ofstream::out);
for(int ii=0;ii<mesh->num_vertices();ii++)
{
    auto actual_vertex=Vertex(*mesh, ii);
 	auto point_vertex=actual_vertex.point(); 
outputFile << point_vertex[0]<<", "<<point_vertex[1]<<", "<<vertex_color[ii]<<", "<<actual_vertex.global_index()<<", "<<ii<<", "<<use_vertex_color[ii]<<", "<<cont<<"\n";
    
}

}
