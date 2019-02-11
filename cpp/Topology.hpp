#include <dolfin.h>
#ifndef TOPOLOGYFENICS_H
#define TOPOLOGYFENICS_H

class Topology{
private:
MeshConnectivity _topology_N2F;
MeshConnectivity _topology_F2N;
MeshConnectivity _topology_K2N;
MeshConnectivity _topology_K2F;
MeshConnectivity _topology_F2K;
MeshConnectivity _topology_N2K;
MeshConnectivity _topology_N2E;
MeshConnectivity _topology_E2N;
MeshConnectivity _topology_N2N;

MeshConnectivity init_topology(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &from,const unsigned int &to);
public:
Topology(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &gdim);// constructor
~Topology(void){ // CHECK THIS
_topology_N2F.clear();
_topology_F2N.clear();
_topology_K2N.clear();
_topology_K2F.clear();
_topology_F2K.clear();
_topology_N2K.clear();
_topology_N2E.clear();
_topology_E2N.clear();
_topology_N2N.clear();
};
MeshConnectivity N2F(void){return _topology_N2F;};
MeshConnectivity F2N(void){return _topology_F2N;};
MeshConnectivity K2N(void){return _topology_K2N;};
MeshConnectivity K2F(void){return _topology_K2F;};
MeshConnectivity F2K(void){return _topology_F2K;};
MeshConnectivity N2K(void){return _topology_N2K;};
MeshConnectivity N2E(void){return _topology_N2E;};
MeshConnectivity E2N(void){return _topology_E2N;};
MeshConnectivity N2N(void){return _topology_N2N;};
MeshConnectivity topologyN2N(std::shared_ptr<dolfin::Mesh> mesh, dolfin::MeshConnectivity &topology_N2F, dolfin::MeshConnectivity &topology_F2N);
std::vector<std::vector< unsigned int > > N2PatchDofs( std::shared_ptr<dolfin::Mesh> mesh,std::shared_ptr<const dolfin::GenericDofMap> dofmap, const dolfin::MeshConnectivity &topology_N2F);

};

MeshConnectivity Topology::init_topology(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &from,const unsigned int &to )
{
mesh->init(from,to);
return mesh->topology()(from,to);}



Topology::Topology(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &gdim)
:
_topology_N2F(init_topology(mesh,0,gdim-1)),
_topology_F2N(init_topology(mesh,gdim-1,0)),
_topology_K2N(init_topology(mesh,gdim,0)),
_topology_K2F(init_topology(mesh,gdim,gdim-1)),
_topology_F2K(init_topology(mesh,gdim-1,gdim)),
_topology_N2K(init_topology(mesh,0,gdim)),
_topology_N2E(init_topology(mesh,0,1)),  // we consider N2E, E2N to compute N2N also for the 3D case 
_topology_E2N(init_topology(mesh,1,0)),  // we consider N2E, E2N to compute N2N also for the 3D case 
_topology_N2N(topologyN2N(mesh,_topology_N2E,_topology_E2N))
{}





MeshConnectivity Topology::topologyN2N( std::shared_ptr<dolfin::Mesh> mesh, dolfin::MeshConnectivity &topology_N2E, dolfin::MeshConnectivity &topology_E2N)
{
    MeshConnectivity _topologyN2N(0,0);
    std::vector< std::vector<std::size_t> > topology_N2N(mesh->num_vertices());
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
	  for (MeshEntityIterator ee(  MeshEntity(*mesh, 0, nn),1 ); !ee.end(); ++ee)
	 {
		auto edge_dof=topology_N2E(nn)[ee.pos()];
		// loop on all the vertices belonging to the edge 
		for(int nn2=0;nn2<gdim;nn2++)
		   // if the node is different from the actual one, add it to its patch
		   if(nn!=topology_E2N(edge_dof)[nn2])
			topology_N2N[nn].push_back(topology_E2N(edge_dof)[nn2]);
           
	 }
	}
	
	_topologyN2N.set(topology_N2N);
	return _topologyN2N;
}   





std::vector<std::vector< unsigned int > > Topology::N2PatchDofs( std::shared_ptr<dolfin::Mesh> mesh,std::shared_ptr<const dolfin::GenericDofMap> dofmap, const dolfin::MeshConnectivity &topology_N2F)
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

#endif