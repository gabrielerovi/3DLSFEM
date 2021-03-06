from dolfin import *
import numpy as np
import GaussPointsWeights as GPW 
from scipy.sparse import lil_matrix


DOLFIN_TOLL=0.000000000001



def AreaTriangle(face_vertices): 
 Area= (0.5 * np.linalg.norm(np.cross(face_vertices[1]-face_vertices[0],face_vertices[2]-face_vertices[0])))
 return Area   

def VolumeTetrahedron(tetrahedron_vertices): 
 Volume=1.0/6.0*np.dot(tetrahedron_vertices[0]-tetrahedron_vertices[3], np.cross(tetrahedron_vertices[1]-tetrahedron_vertices[3],tetrahedron_vertices[2]-tetrahedron_vertices[3]))
 Volume=np.absolute(Volume)
 return Volume  
 
 
 
# here normal_direction=1 if the normal of the face is cross(FV(1)-FV(0),FV(2)-FV(0))
# here normal_direction=-1 if the normal of the face is cross(-FV(1)+FV(0),-FV(2)+FV(0))

def RaviartThomas3D(face_vertices,opposite_vertex,q_points,normal_direction):

 dim=3
 qp_length=np.size(q_points,0)
 RT=np.zeros((qp_length,dim))
 midpoint=np.mean(face_vertices, axis=0) 
 tmp=midpoint-opposite_vertex
 normal_face=np.cross( face_vertices[1]-face_vertices[0] ,face_vertices[2]-face_vertices[0] )
 normal_face=normal_face*normal_direction

 shape_function_sign=np.sign(np.dot(normal_face,tmp))
 tetrahedron_vertices=np.zeros((4,3))
 tetrahedron_vertices[0:3]=face_vertices
 tetrahedron_vertices[3]=opposite_vertex

 RT_coeff=shape_function_sign/(dim*VolumeTetrahedron(tetrahedron_vertices))
 # loop on all the quadrature points
 for qp in range(0,qp_length):
    for dd in range(0,dim):
        RT[qp,dd]=RT_coeff*(q_points[qp,dd]-opposite_vertex[dd])
        
 return RT



def fine_face_integral_of_coarse_function(fine_face_vertices,normal_fine,coarse_face_vertices,normal_coarse,coarse_opposite_vertex,normal_direction ):
  qrule=4
  dim=3
  q_points_ref,w_qp=GPW.TriGaussPoints2D(qrule)   
  qp_length=np.size(q_points_ref,0)   
  integral = 0.0
  J=np.array([ [fine_face_vertices[1,0]-fine_face_vertices[0,0], fine_face_vertices[2,0]-fine_face_vertices[0,0]],
               [fine_face_vertices[1,1]-fine_face_vertices[0,1], fine_face_vertices[2,1]-fine_face_vertices[0,1]],
               [fine_face_vertices[1,2]-fine_face_vertices[0,2], fine_face_vertices[2,2]-fine_face_vertices[0,2]] ])
  bJ=np.array([fine_face_vertices[0,0],fine_face_vertices[0,1],fine_face_vertices[0,2]])
  q_points=np.zeros((qp_length,dim))
  for qp in range(0,qp_length):
      q_points[qp,:]=bJ+J.dot(q_points_ref[qp])

  RT3D=RaviartThomas3D(coarse_face_vertices,coarse_opposite_vertex,q_points,normal_direction)
  for qp in range(0,qp_length):
   integral= integral + w_qp[qp]*np.dot(RT3D[qp],normal_fine)
  integral = integral * AreaTriangle(fine_face_vertices) 
  return integral  
  
  
  
  
  
def mixed_interpolation_matrix(mesh_coarse, mesh_fine, WC,WF, dofmap_coarse, dofmap_fine, mesh_topologyC, mesh_topologyF):
  # consider a face with vertices n1<n2<n3
  # if normal_direction=1:  cross(n2-n1,n3-n1) outward
  # if normal_direction=-1:  cross(n2-n1,n3-n1) inward
  normal_direction=-1
  gdim=mesh_coarse.geometry().dim()
  num_of_subdivision=2**(gdim)
  WC_num_dofs=WC.dim()
  WF_num_dofs=WF.dim();
  WC2WFnumpy=lil_matrix((WF_num_dofs, WC_num_dofs))

  mesh_points_coarse=mesh_coarse.coordinates()
  mesh_points_fine=mesh_fine.coordinates()
  
  N_2_F=mesh_topologyC.N_2_F 
  F_2_N=mesh_topologyC.F_2_N
  K_2_N=mesh_topologyC.K_2_N
  K_2_F=mesh_topologyC.K_2_F
    
  N_2_F_fine=mesh_topologyF.N_2_F
  F_2_N_fine=mesh_topologyF.F_2_N
  K_2_N_fine=mesh_topologyF.K_2_N
  K_2_F_fine=mesh_topologyF.K_2_F   
    
  # build a local vector which, given the faces, returns the opposite node
  coarse_node_opposite_to_face=np.zeros(gdim+1, dtype=int)
  # loop on all the coarse cells
  # we denote by K the coarse element
  for cell in cells(mesh_coarse):
      # take the coarse cell index
      cell_index=cell.index() 
      # F_of_K_fine, F_of_K_coarse: coarse and fine faces
      F_of_K_fine=np.empty(0, dtype=int)
      F_of_K_coarse= K_2_F(cell.index())
      # N_of_K_fine, N_of_K_coarse: coarse and fine nodes      
      N_of_K_fine=np.empty(0, dtype=int)
      N_of_K_coarse= K_2_N(cell.index())   
      # loop on all the fine elements inside the coarse element cell.index()
      # find all the fine nodes and faces
      for fine in range(num_of_subdivision*cell_index,num_of_subdivision*(cell_index+1)):
         N_of_K_fine=np.concatenate((N_of_K_fine,K_2_N_fine(fine)),axis=0)
         F_of_K_fine=np.concatenate((F_of_K_fine,K_2_F_fine(fine)),axis=0)
      # some nodes/edges are shared among the fine elements, therefore remove the duplicates   
      N_of_K_fine=np.unique( N_of_K_fine )  
      F_of_K_fine=np.unique( F_of_K_fine )  
      
      # consider only the nodes which are on the fine mesh, but not on the coarse one
      only_fine_nodes=np.setdiff1d(N_of_K_fine, N_of_K_coarse)
      
      # loop on coarse faces to find, given the coarse face ff1, the opposite coarse vertex
      for ff1 in range(0,F_of_K_coarse.size):
          # compute the face nodes of the coarse face
          coarse_face_nodes=np.sort(F_2_N(F_of_K_coarse[ff1]))
          # loop on the coarse nodes and see which one is not in coarse_face_nodes: that one is the opposite node
          for nn1 in range(0,N_of_K_coarse.size):
           counter=0
           for dd in range(0,gdim):
            if(N_of_K_coarse[nn1]!=coarse_face_nodes[dd]):
              counter=counter+1
           if(counter==3):
              coarse_node_opposite_to_face[ff1]=N_of_K_coarse[nn1]
              break
        
          
      # define the coordinates relative to N_of_K_coarse
      coordinates_of_K_coarse=mesh_points_fine[N_of_K_coarse]
      # define the coordinates relative to N_of_K_only_fine 
      coordinates_of_K_only_fine=mesh_points_fine[np.setdiff1d(N_of_K_fine, N_of_K_coarse)]
      
      # compute the coarse midpoint of the coarse element 
      
      # CERCA PER NUMPY MEAN O SIMILI
      # CONTROLLA SE NODO 0 E' OPPOSTO A FACCIA N1-N2-N3
      mean_coordinates_of_K_coarse=np.matrix(
      [(coordinates_of_K_coarse[1]+coordinates_of_K_coarse[2]+coordinates_of_K_coarse[3] ) /3.0,
      (coordinates_of_K_coarse[0]+coordinates_of_K_coarse[2]+coordinates_of_K_coarse[3] ) /3.0,  
      (coordinates_of_K_coarse[0]+coordinates_of_K_coarse[1]+coordinates_of_K_coarse[3] ) /3.0,
      (coordinates_of_K_coarse[0]+coordinates_of_K_coarse[1]+coordinates_of_K_coarse[2] ) /3.0       ])
      
      mean_edge_coordinates_of_K_coarse=np.array([
      [(coordinates_of_K_coarse[0]+coordinates_of_K_coarse[1])*0.5],
      [(coordinates_of_K_coarse[0]+coordinates_of_K_coarse[2])*0.5],
      [(coordinates_of_K_coarse[0]+coordinates_of_K_coarse[3])*0.5],
      [(coordinates_of_K_coarse[1]+coordinates_of_K_coarse[2])*0.5],
      [(coordinates_of_K_coarse[1]+coordinates_of_K_coarse[3])*0.5],
      [(coordinates_of_K_coarse[2]+coordinates_of_K_coarse[3])*0.5]
      ])
      mean_edge_indices_of_K_coarse=np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
      # consider local vertices of a given edge
      edge_vertices=np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]])
      # dofs of the coarse nodes
      coarse_node_dofs=dofmap_coarse.entity_dofs(mesh_coarse,0,np.array(N_of_K_coarse,dtype=np.uintp))
      # dofs of the only fine nodes
      only_fine_dofs=dofmap_fine.entity_dofs(mesh_fine,0,np.array(only_fine_nodes,dtype=np.uintp))
      # dofs on all the nodes
      coarse_and_fine_dofs=dofmap_fine.entity_dofs(mesh_fine,0,np.array(N_of_K_coarse,dtype=np.uintp))
      # loop on the only-fine vertices
      for nn1 in range(0,only_fine_nodes.size):
           # loop on the edges:
           for ee2 in range(0,(gdim-1)*3):
            if(np.linalg.norm(np.subtract(coordinates_of_K_only_fine[nn1],mean_edge_coordinates_of_K_coarse[ee2]))<DOLFIN_TOLL):
                  tmp=dofmap_fine.entity_dofs(mesh_fine,0,np.array([only_fine_nodes[nn1]],dtype=np.uintp))
                  edge_vertices=N_of_K_coarse[mean_edge_indices_of_K_coarse[ee2]]
                  tmpC1=dofmap_coarse.entity_dofs(mesh_coarse,0,np.array([edge_vertices[0]],dtype=np.uintp))
                  tmpC2=dofmap_coarse.entity_dofs(mesh_coarse,0,np.array([edge_vertices[1]],dtype=np.uintp))
                  for dd in range(0,gdim):
                      WC2WFnumpy[tmp[dd],[tmpC1[dd], tmpC2[dd]]]=[0.5, 0.5]
                  break
 
            
           
#            loop on the coarse vertices
#            for nn2 in range(0,N_of_K_coarse.size):
#                     check if the fine vertex is actually the midpoint of two coarse vertices
#                     if(np.linalg.norm(np.subtract(coordinates_of_K_only_fine[nn1],mean_coordinates_of_K_coarse[nn2]))<DOLFIN_TOLL):
#                         tmp=dofmap_fine.entity_dofs(mesh_fine,0,np.array([only_fine_nodes[nn1]],dtype=np.uintp))
#                         tmpC1=dofmap_coarse.entity_dofs(mesh_coarse,0,np.array([N_of_K_coarse[edge_vertices.item((nn2,0))]],dtype=np.uintp))
#                         tmpC2=dofmap_coarse.entity_dofs(mesh_coarse,0,np.array([N_of_K_coarse[edge_vertices.item((nn2,1))]],dtype=np.uintp))
#                         WC2WFnumpy[tmp[0],[tmpC1[0], tmpC2[0]]]=[0.5, 0.5]
#                         WC2WFnumpy[tmp[1],[tmpC1[1], tmpC2[1]]]=[0.5, 0.5]
      # the dofs related to the nodes which are on both the coarse and fine meshes remain the same
      
      # CONTROLLA QUI. GIUSTO CHE SIANO COARSE AND FINE? NON SOLO COARSE?
      # magari il nome e' misleading, perche' cmq sono i nodi coarse ma enumerati nella mesh fine 

      for nn1 in range(0,  coarse_and_fine_dofs.size):
       WC2WFnumpy[coarse_and_fine_dofs[nn1],coarse_node_dofs[nn1]]=1
       
      # now we consider the interpolation for the RAVIART-THOMAS
      # loop on coarse edges
      for ff1 in range(0,F_of_K_coarse.size): 
          # given a coarse face ff1, consider its gdim coarse nodes: sort them, to have a unique representation of the normal (as in Marie Rognes)
          coarse_face_nodes=np.sort(F_2_N(F_of_K_coarse[ff1]))
          # compute the nodes of the face, with n1 < n2 < n3
          coarse_face_vertices=mesh_points_coarse[coarse_face_nodes]
          # compute the normal of the face 
          # do not need to normalize because then we just consider the sign
          normal_ff1=np.cross( coarse_face_vertices[1]-coarse_face_vertices[0] ,coarse_face_vertices[2]-coarse_face_vertices[0] ) 
          normal_ff1=normal_ff1*normal_direction
          area_ff1=np.linalg.norm(normal_ff1*0.5)
          normal_ff1=normal_ff1/np.linalg.norm(normal_ff1)
          #normal_ff1=normal_ff1/np.linalg.norm(normal_ff1)
          # consider the coarse face and compute the vector that goes from the opposite coarse node to the midpoint of the face
          face_coarse=Facet(mesh_coarse,F_of_K_coarse[ff1])
          coarse_face_midpoint=np.array([face_coarse.midpoint().x(),face_coarse.midpoint().y(),face_coarse.midpoint().z()])
          midpoint_minus_opposite_node=coarse_face_midpoint-mesh_points_coarse[coarse_node_opposite_to_face[ff1]]

          # compute the sign of the coarse function
#           print("sign_of_coarse_function")
#           print(coarse_node_opposite_to_face)
#           print(coarse_face_midpoint)
#           print(coarse_face_vertices)
#           print(mesh_points_coarse[coarse_node_opposite_to_face[ff1]])
#           print()
#           print(midpoint_minus_opposite_node,normal_ff1,np.dot(midpoint_minus_opposite_node,normal_ff1))
          sign_of_coarse_function=np.sign(np.dot(midpoint_minus_opposite_node,normal_ff1))
          # loop on fine edges 
          for ff2 in range(0,F_of_K_fine.size):    
           face=Facet(mesh_fine,F_of_K_fine[ff2])
           midpoint_fine=np.array([face.midpoint().x(),face.midpoint().y()])
           # fine nodes of the given fine edge ff2
           fine_face_nodes=np.sort(F_2_N_fine(F_of_K_fine[ff2]))
           # coordinates of the fine nodes
           fine_face_vertices=mesh_points_fine[fine_face_nodes]
           # fine normal
           
           ############## QUI USI IL FATTO CHE SIA 3D
           normal_ff2=np.cross( fine_face_vertices[1]-fine_face_vertices[0] ,fine_face_vertices[2]-fine_face_vertices[0] )
           normal_ff2=normal_ff2*normal_direction
           normal_ff2=normal_ff2/np.linalg.norm(normal_ff2)
           ############## FAI IL CROSS PRODUCT PER BBENE ###################
           fine_face_midpoint=np.array([Facet(mesh_fine,F_of_K_fine[ff2]).midpoint().x(),Facet(mesh_fine,F_of_K_fine[ff2]).midpoint().y(),Facet(mesh_fine,F_of_K_fine[ff2]).midpoint().z()])
           tmp=fine_face_midpoint-mesh_points_coarse[coarse_node_opposite_to_face[ff1]]
           # compute the dot product between the fine normal and the vector that goes from the opposite coarse node to the midpoint of the fine face
           dot_product=np.dot(tmp,normal_ff2)
           # if it is very small, then put it to zero (they are normal)                      
           if(np.absolute(dot_product)<DOLFIN_TOLL):
              sign_dot_product=0.0
           else:
           # otherwise consider the sign
              sign_dot_product=np.sign(dot_product)
           
           
           # Now we check wether the fine face belongs to the coarse face
           # if it is inside, then value=0.25, otherwise value=0.125
#           if(ff2==0):
#              print("000000000000")
#              print(fine_face_vertices)
#              print(sign_dot_product,sign_of_coarse_function )
           
           
           value=fine_face_integral_of_coarse_function(fine_face_vertices,normal_ff2,coarse_face_vertices,normal_ff1,mesh_points_coarse[coarse_node_opposite_to_face[ff1]],normal_direction)
           dof_face_fine=dofmap_fine.entity_dofs(mesh_fine,gdim-1,np.array([F_of_K_fine[ff2]],dtype=np.uintp))
           dof_face_coarse=dofmap_coarse.entity_dofs(mesh_coarse,gdim-1,np.array([F_of_K_coarse[ff1]],dtype=np.uintp))
#            print("VALUE",value)
#            print("fine_face_vertices",fine_face_vertices)
#            print("normal_ff2",normal_ff2)
#            print("coarse_face_vertices",coarse_face_vertices)
#            print("ff1",ff1)
#            print("ff2",ff2)
#            print("normal_ff1",normal_ff1)
#            print("mesh_points_coarse[coarse_node_opposite_to_face[ff1]]",mesh_points_coarse[coarse_node_opposite_to_face[ff1]])
#            print()
           
#            prova=np.array([[0,0,0.5],[0.5,0.5,0],[0,0.5,0.5]])
#            prova_coarse=np.array([[1,0,.0],[0.,1.5,0],[0,0.,2]])
#            opposite_prova=np.array([0,0,0])
#            normal_prova=np.array([0.7071,0,0.7071])
#            normal_prova=normal_prova/np.linalg.norm(normal_prova)
#            normal_prova_coarse= np.array([0.5774,    0.5774,    0.5774]) 
#            q_points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
#            ok=RT.RaviartThomas3D(prova,opposite_prova,q_points,normal_direction)
#            value=fine_face_integral_of_coarse_function(prova,normal_prova,prova_coarse,normal_prova_coarse,opposite_prova,normal_direction)
#            print("ora")
#            print(value)
#            print()
#            print()
#            print(value[1])
           for dd in range(0,gdim):
               WC2WFnumpy[dof_face_fine[dd],dof_face_coarse[dd]]=value
           value=0.25
#            for mm1 in range(0,gdim):
#                tmp=fine_face_vertices[mm1]-coarse_face_vertices[0]
#                # if at least one fine point is not on the plane, exit
# #                if(ff2==0): 
# #                   print("np.absolute(np.dot(normal_ff1,tmp))")
# #                   print(np.absolute(np.dot(normal_ff1,tmp)))
#                if(np.absolute(np.dot(normal_ff1,tmp)) >DOLFIN_TOLL):
#                   value=0.125
#                   break
#                # otherwise check if it is inside the triangle
#                else:
#                   # if it is not, break; otherwise, continue
#                   alpha=np.linalg.norm(np.cross(fine_face_vertices[mm1] - coarse_face_vertices[1], fine_face_vertices[mm1] -coarse_face_vertices[2]))/(2*area_ff1)
#                   beta=np.linalg.norm(np.cross(fine_face_vertices[mm1] - coarse_face_vertices[2] , fine_face_vertices[mm1] -coarse_face_vertices[0]))/(2*area_ff1)
#                   gamma=1-beta-alpha
# #                   if(ff2==0): 
# #                     print("alpha,beta,gamma")
# #                     print(alpha,beta,gamma)
#                   if(alpha<=1+DOLFIN_TOLL and beta <=1+DOLFIN_TOLL and -DOLFIN_TOLL<=gamma and gamma<=1+DOLFIN_TOLL):
#                        value=0.25
#                   else:
#                        value=0.125
#                        break
# 
#            # the degrees of freedom
#            dof_face_fine=dofmap_fine.entity_dofs(mesh_fine,gdim-1,np.array([F_of_K_fine[ff2]],dtype=np.uintp))
#            dof_face_coarse=dofmap_coarse.entity_dofs(mesh_coarse,gdim-1,np.array([F_of_K_coarse[ff1]],dtype=np.uintp))
#            for dd in range(0,gdim):
#               WC2WFnumpy[dof_face_fine[dd],dof_face_coarse[dd]]=sign_dot_product*value*sign_of_coarse_function  
#            if(ff2==11):
#             print("ff1",ff1)
#             print("value")              
#             print(value,sign_dot_product*value*sign_of_coarse_function )
#            print(sign_dot_product)
#            print(sign_of_coarse_function)
#            print(WC2WFnumpy[dof_face_fine[dd],dof_face_coarse[dd]])


  return WC2WFnumpy

