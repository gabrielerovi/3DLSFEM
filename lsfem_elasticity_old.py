from dolfin import *
from scipy.io import savemat
from petsc4py import PETSc
from scipy.sparse import csr_matrix
import numpy as np

def boundary(x, on_boundary):
    return on_boundary
def A(sigma):
	return beta*sigma+alpha*tr(sigma)*Identity(2)
def epsilon(u):
    return sym(grad(u))
    
##########################################################
################ COARSE MESH #############################  
##########################################################  
aL = 0.0
aR = 1.0
bB = 0.0
bT = 1.0
Na = 3
Nb = 3
mesh = RectangleMesh(Point(aL, bB),Point(aR, bT), Na,Nb) #,"crossed"
gdim = mesh.geometry().dim()
mesh_points=mesh.coordinates()
N=mesh.num_vertices()

######################################################################
################ COARSE MESH TOPOLOGY INFORMATION #################### 
######################################################################   
# node to edge        
mesh.init(0,1)
N_2_F=mesh.topology()(0,1)
# edge to node     
mesh.init(1,0)
F_2_N=mesh.topology()(1,0)
# element to node
mesh.init(gdim,0)
K_2_N=mesh.topology()(gdim,0)
# element to edge
mesh.init(gdim,1)
K_2_E=mesh.topology()(gdim,1)

##########################################################
################ FINE  MESH  #############################  
##########################################################  
mesh2=refine(mesh)
mesh_points2=mesh2.coordinates()

######################################################################
################ FINE  MESH TOPOLOGY  INFORMATION #################### 
######################################################################  
# node to edge fine       
mesh2.init(0,1)
N_2_F_fine=mesh2.topology()(0,1)
# edge to node fine  
mesh2.init(1,0)
F_2_N_fine=mesh2.topology()(1,0)
# element to node fine
mesh2.init(gdim,0)
K_2_N_fine=mesh2.topology()(gdim,0)
# element to edge fine
mesh2.init(gdim,1)
K_2_E_fine=mesh2.topology()(gdim,1)

#################################################    
################ PARAMETERS #####################  
#################################################  

# Functioncal constants: F=C_equilibrium||div sigma +f ||^2  + C_constitutive ||A(sigma)-epsilon(u) ||^2
C_equilibrium=1
C_constitutive=1.0
# Lame parameters and other related constants
Lambda = Constant(1.0)
Mu = Constant(1.0)
beta=Constant(1.0/(2*Mu))
alpha = Constant(-beta * Lambda /(2*Lambda +2*Mu))
gamma=Constant(Lambda+2*Mu);
# used tollerance
DOLFIN_TOLL=0.000000000001

# external force f and dirichlet boundary displacement g
f = Expression(("( 4.0 * pi*pi * sin(pi * x[0]) * sin(pi * x[1]) - 2.0 * pi*pi *cos( pi * x[0]) * cos( pi * x[1]) )","( 4.0 * pi*pi * sin(pi * x[0]) * sin(pi * x[1]) - 2.0 * pi*pi *cos( pi * x[0]) * cos( pi * x[1]) )"),degree=4)
g = Expression(("0.0","0.0"),degree=4)

###########################################################        
################ COARSE FUNCTION SPACE ####################
###########################################################   
# 2D, 2 components Lagrange FEM  
fe2DP1 = VectorElement("P",mesh.ufl_cell(),1)
# 2D, 2 components RT FEM  
fe2DRT = VectorElement("RT",mesh.ufl_cell(),1)
# 2D, 4 components mixed FEM 
feW  = MixedElement([fe2DRT,fe2DP1])
W = FunctionSpace(mesh, feW) 
w = Function(W)
sigma,u = split(w)
W_num_dofs=W.dim();
dm = W.dofmap()

 
#####################################################################         
################ COARSE PATCH: FROM NODE TO DOFS  ###################
##################################################################### 
Patch_coarse=np.empty((N,),dtype=object)  

# loop on all the nodes
for nn in range(0,N):
   # for the given node, find all the nodal dofs related to it
   tmp=dm.entity_dofs(mesh,0,np.array([nn],dtype=np.uintp))
   # then consider all the edges/faces (2D/3D) connecting the node nn
   N_to_F=N_2_F(nn);
   # for each of these, consider all the dofs inside the face
   for ff1 in range(0,N_to_F.size):
      tmp3=dm.entity_dofs(mesh,1,np.array([N_to_F[ff1]],dtype=np.uintp)) 
      tmp=np.concatenate((tmp,tmp3),axis=0)
   # collect all the nodes of the coarse patch
   tmp=np.unique(tmp)
   Patch_coarse[nn]=tmp





###########################################################        
################ FINE  FUNCTION  SPACE ####################
########################################################### 
# 2D, 2 components Lagrange FEM  
fe2DP12 = VectorElement("P",mesh2.ufl_cell(),1)
# 2D, 2 components RT FEM  
fe2DRT2 = VectorElement("RT",mesh2.ufl_cell(),1)
# 2D, 4 components mixed FEM 
feW2  = MixedElement([fe2DRT2,fe2DP12])
W2 = FunctionSpace(mesh2, feW2) #MixedFunctionSpace([RT,P1])
dm2 = W2.dofmap()
W2_num_dofs=W2.dim();




#################################################################################                
################ INTERPOLATION OPERATOR MIXED FUNCTION SPACE ####################
################################################################################# 
# Create the interpolation operator from WC (Coarse) to WF (Fine)
WC2WFnumpy=csr_matrix((W2_num_dofs, W_num_dofs))

# build a local vector which, given the edges, returns the opposite node
coarse_node_opposite_to_edge=np.zeros(gdim+1, dtype=int)

# loop on all the coarse cells
# we denote by K the coarse element
for cell in cells(mesh):
      # take the coarse cell index
      cell_index=cell.index() 
      # E_of_K_fine, E_of_K_coarse: coarse and fine edges
      E_of_K_fine=np.empty(0, dtype=int)
      E_of_K_coarse= K_2_E(cell.index())
      # N_of_K_fine, N_of_K_coarse: coarse and fine nodes      
      N_of_K_fine=np.empty(0, dtype=int)
      N_of_K_coarse= K_2_N(cell.index())   
      # loop on all the fine elements inside the coarse element cell.index()
      # find all the fine nodes and edges
      for fine in range(4*cell_index,4*(cell_index+1)):
         N_of_K_fine=np.concatenate((N_of_K_fine,K_2_N_fine(fine)),axis=0)
         E_of_K_fine=np.concatenate((E_of_K_fine,K_2_E_fine(fine)),axis=0)
      # some nodes/edges are shared among the fine elements, therefore remove the duplicates   
      N_of_K_fine=np.unique( N_of_K_fine )  
      E_of_K_fine=np.unique( E_of_K_fine )  
      
      # consider only the nodes which are on the fine mesh, but not on the coarse one
      only_fine_nodes=np.setdiff1d(N_of_K_fine, N_of_K_coarse)
      
      # loop on coarse edges to find, given the coarse edge ee1, the opposite coarse vertex
      for ee1 in range(0,E_of_K_coarse.size):
          # compute the extreme two nodes of the coarse edge
          coarse_edge_nodes=np.sort(F_2_N(E_of_K_coarse[ee1]))
          # loop on the coarse nodes and see which one is not in coarse_edge_nodes: that one is the opposite node
          for nn1 in range(0,N_of_K_coarse.size):
           if(N_of_K_coarse[nn1]!=coarse_edge_nodes[0] and N_of_K_coarse[nn1]!=coarse_edge_nodes[1]):
             coarse_node_opposite_to_edge[ee1]=N_of_K_coarse[nn1]
                        
      # define the coordinates relative to N_of_K_coarse
      coordinates_of_K_coarse=mesh_points2[N_of_K_coarse]
      # define the coordinates relative to N_of_K_only_fine 
      coordinates_of_K_only_fine=mesh_points2[np.setdiff1d(N_of_K_fine, N_of_K_coarse)]
      
      # compute the coarse midpoint of the coarse element 
      mean_coordinates_of_K_coarse=np.matrix(
      [(coordinates_of_K_coarse[1]+coordinates_of_K_coarse[2] ) /2,
      (coordinates_of_K_coarse[0]+coordinates_of_K_coarse[2] ) /2,  
      (coordinates_of_K_coarse[0]+coordinates_of_K_coarse[1] ) /2 ])
      # consider local vertices of a given edge
      edge_vertices=np.matrix([[1,2],[0,2],[0,1]])
      # dofs of the coarse nodes
      coarse_node_dofs=dm.entity_dofs(mesh,0,np.array(N_of_K_coarse,dtype=np.uintp))
      # dofs of the only fine nodes
      only_fine_dofs=dm2.entity_dofs(mesh2,0,np.array(only_fine_nodes,dtype=np.uintp))
      # dofs on all the nodes
      coarse_and_fine_dofs=dm2.entity_dofs(mesh2,0,np.array(N_of_K_coarse,dtype=np.uintp))
      # loop on the only-fine vertices
      for nn1 in range(0,only_fine_nodes.size):
           # loop on the coarse vertices
           for nn2 in range(0,N_of_K_coarse.size):
                    #check if the fine vertex is actually the midpoint of two coarse vertices
                    if(np.linalg.norm(np.subtract(coordinates_of_K_only_fine[nn1],mean_coordinates_of_K_coarse[nn2]))<DOLFIN_TOLL):
                        tmp=dm2.entity_dofs(mesh2,0,np.array([only_fine_nodes[nn1]],dtype=np.uintp))
                        tmpC1=dm.entity_dofs(mesh,0,np.array([N_of_K_coarse[edge_vertices.item((nn2,0))]],dtype=np.uintp))
                        tmpC2=dm.entity_dofs(mesh,0,np.array([N_of_K_coarse[edge_vertices.item((nn2,1))]],dtype=np.uintp))
                        WC2WFnumpy[tmp[0],[tmpC1[0], tmpC2[0]]]=[0.5, 0.5]
                        WC2WFnumpy[tmp[1],[tmpC1[1], tmpC2[1]]]=[0.5, 0.5]
      # the dofs related to the nodes which are on both the coarse and fine meshes remain the same
      for nn1 in range(0,  coarse_and_fine_dofs.size):
       WC2WFnumpy[coarse_and_fine_dofs[nn1],coarse_node_dofs[nn1]]=1
       
      # now we consider the interpolation for the raviart-thomas 
      # loop on coarse edges
      for ee1 in range(0,E_of_K_coarse.size): 
          # given a coarse edge ee1, consider its two coarse nodes: sort them, to have a unique representation of the normal  (as in Marie Rognes)
          coarse_edge_nodes=np.sort(F_2_N(E_of_K_coarse[ee1]))
          # compute the nodes of the edge, with n1 < n2, being tangent vector that goes from n1 to n2
          coarse_edge_vertices=mesh_points[coarse_edge_nodes]
          # compute the normal of the edge (as in Marie Rognes, is the clockwise rotation of the tangent vector from n1 to n2)
          # do not need to normalize because then we just consider the sign
          normal_ee1=np.array([coarse_edge_vertices.item(3)-coarse_edge_vertices.item(1),coarse_edge_vertices.item(0)-coarse_edge_vertices.item(2) ]) 
          #normal_ee1=normal_ee1/np.linalg.norm(normal_ee1)
          # consider the coarse face of the given edge and compute the vector that goes from the opposite coarse node to the midpoint of the face
          face_coarse=Facet(mesh,E_of_K_coarse[ee1])
          midpoint_minus_opposite_node=np.array([face_coarse.midpoint().x()-mesh_points[coarse_node_opposite_to_edge[ee1]][0],face_coarse.midpoint().y()-mesh_points[coarse_node_opposite_to_edge[ee1]][1]])
          # compute the sign of the coarse function
          sign_of_coarse_function=np.sign(np.dot(midpoint_minus_opposite_node,normal_ee1))
          # loop on fine edges 
          for ee2 in range(0,E_of_K_fine.size):    
           face=Facet(mesh2,E_of_K_fine[ee2])
           midpoint_fine=np.array([face.midpoint().x(),face.midpoint().y()])
           # fine nodes of the given fine edge ee2
           fine_edge_nodes=np.sort(F_2_N_fine(E_of_K_fine[ee2]))
           # coordinates of the fine nodes
           fine_edge_vertices=mesh_points2[fine_edge_nodes]
           # fine normal
           normal_ee2=np.array([fine_edge_vertices.item(3)-fine_edge_vertices.item(1),fine_edge_vertices.item(0)-fine_edge_vertices.item(2) ]) 
           fine_edge_midpoint=np.array([Facet(mesh2,E_of_K_fine[ee2]).midpoint().x(),Facet(mesh2,E_of_K_fine[ee2]).midpoint().y()])
           tmp=fine_edge_midpoint-mesh_points[coarse_node_opposite_to_edge[ee1]]
           # compute the dot product between the fine normal and the vector that goes from the opposite coarse node to the midpoint of the fine edge
           dot_product=np.dot(tmp,normal_ee2)
           # if it is very small, then put it to zero (they are normal)                      
           if(np.absolute(dot_product)<DOLFIN_TOLL):
              sign_dot_product=0.0
           else:
           # otherwise consider the sign
              sign_dot_product=np.sign(dot_product)
           
           # We have two different cases of absolute values of the interpolation operator, i.e. a) and b) 
           # a) all the fine edges which are on the boundary of the coarse element and on the median
           # b) all the fine edges whose extreme nodes do not belong to the nodes of the coarse mesh
           # if all try_i=0, then we are in the a) case
           # otherwise in the b)       
           # Node_E_C_1 - Node_E_F_1
           try1=np.linalg.norm(np.array([coarse_edge_vertices.item(0)-fine_edge_vertices.item(0),coarse_edge_vertices.item(1)-fine_edge_vertices.item(1) ]))
           # Node_E_C_1 - Node_E_F_2
           try2=np.linalg.norm(np.array([coarse_edge_vertices.item(0)-fine_edge_vertices.item(2),coarse_edge_vertices.item(1)-fine_edge_vertices.item(3) ]))           
           # Node_E_C_2 - Node_E_F_1
           try3=np.linalg.norm(np.array([coarse_edge_vertices.item(2)-fine_edge_vertices.item(0),coarse_edge_vertices.item(3)-fine_edge_vertices.item(1) ]))
           # Node_E_C_2 - Node_E_F_2
           try4=np.linalg.norm(np.array([coarse_edge_vertices.item(2)-fine_edge_vertices.item(2),coarse_edge_vertices.item(3)-fine_edge_vertices.item(3) ]))
           # the degrees of freedom
           efine=dm2.entity_dofs(mesh2,1,np.array([E_of_K_fine[ee2]],dtype=np.uintp))
           ecoarse=dm.entity_dofs(mesh,1,np.array([E_of_K_coarse[ee1]],dtype=np.uintp))
           # we change the actual value 
           #if(np.linalg.norm(WC2WFnumpy[efine[0],ecoarse[0]])<DOLFIN_TOLL):
             # a) case with 0.5
           if(try1<DOLFIN_TOLL or try2 < DOLFIN_TOLL or try3 < DOLFIN_TOLL or try4 < DOLFIN_TOLL):
              value=sign_dot_product*0.5*sign_of_coarse_function          
              WC2WFnumpy[efine[0],ecoarse[0]]=value
              WC2WFnumpy[efine[1],ecoarse[1]]=value
             # b) case with 0.25
           else:
              value=sign_dot_product*0.25*sign_of_coarse_function
              WC2WFnumpy[efine[0],ecoarse[0]]=value
              WC2WFnumpy[efine[1],ecoarse[1]]=value
  

    

    

    
J =C_constitutive * inner(A(sigma)-epsilon(u),A(sigma)-epsilon(u))*dx + C_equilibrium * inner(div(sigma)+f,div(sigma)+f)*dx
G = derivative(J,w)


bcs=DirichletBC(W.sub(1), g, DomainBoundary())


      
      
solve(G==0,w,bcs)





sigmacoarse,ucoarse = w.split()



ufile_pvd = File("lsfem_elasticity_u_coarse.pvd")
ufile_pvd << ucoarse
ufile_pvd = File("lsfem_elasticity_sigma_coarse.pvd")
ufile_pvd << sigmacoarse

solutioncoarse=w.vector().array()
solutioncoarse=solutioncoarse[:,None]
#WC2WFnumpy=csr_matrix.todense(WC2WFnumpy)
solutionfine=WC2WFnumpy*solutioncoarse
w2 = Function(W2)
w2.vector()[:] = solutionfine
sigma,u = w2.split()
Gmat = assemble(G)
ufile_pvd << u


ufile_pvd = File("lsfem_elasticity_u.pvd")
ufile_pvd << u
ufile_pvd = File("lsfem_elasticity_sigma.pvd")
ufile_pvd << sigma






print( norm(ucoarse, 'H1', mesh))
print( norm(u, 'H1', mesh))
print( norm(sigmacoarse, 'Hdiv', mesh))
print( norm(sigma, 'Hdiv', mesh))




systemmatrix = derivative(G,w)
systemmatrix=assemble(systemmatrix)
#print("systemmatrix.array()")
#print(systemmatrix.array())

savemat("solutioncoarse.mat",{"solc":solutioncoarse})
savemat("solutionfine.mat",{"solf":solutionfine})
savemat("WC2WFnumpy.mat",{"WC2WFnumpy":WC2WFnumpy})
savemat("coordinates_coarse.mat",{"coordinates_coarse":mesh_points})
savemat("mesh_points2.mat",{"mesh_points2":mesh_points2})
coordinates_edge_coarse=np.zeros((mesh.num_facets(),2))
coordinates_edge_fine=np.zeros((mesh2.num_facets(),2))
dofmap_edge_coarse=np.zeros((mesh.num_facets(),2))
dofmap_edge_fine=np.zeros((mesh2.num_facets(),2))
dofmap_node_coarse=np.zeros((mesh.num_vertices(),2))
dofmap_node_fine=np.zeros((mesh2.num_vertices(),2))
for ee in range(0,mesh.num_facets()):
    coordinates_edge_coarse[[ee],[0]]=Facet(mesh,ee).midpoint().x()
    coordinates_edge_coarse[[ee],[1]]=Facet(mesh,ee).midpoint().y()
    dofmap_edge_coarse[[ee],[0,1]]=dm.entity_dofs(mesh,1,np.array([ee],dtype=np.uintp)) 
for ee in range(0,mesh2.num_facets()):
    coordinates_edge_fine[[ee],[0]]=Facet(mesh2,ee).midpoint().x()
    coordinates_edge_fine[[ee],[1]]=Facet(mesh2,ee).midpoint().y()
    dofmap_edge_fine[[ee],[0,1]]=dm2.entity_dofs(mesh2,1,np.array([ee],dtype=np.uintp)) 
    
for nn in range(0,mesh.num_vertices()):
    dofmap_node_coarse[[nn],[0,1]]=dm.entity_dofs(mesh,0,np.array([nn],dtype=np.uintp))

for nn in range(0,mesh2.num_vertices()):
    dofmap_node_fine[[nn],[0,1]]=dm2.entity_dofs(mesh2,0,np.array([nn],dtype=np.uintp))
savemat("coordinates_edge_coarse.mat",{"coordinates_edge_coarse":coordinates_edge_coarse})
savemat("coordinates_edge_fine.mat",{"coordinates_edge_fine":coordinates_edge_fine})
savemat("dofmap_edge_coarse.mat",{"dofmap_edge_coarse":dofmap_edge_coarse})
savemat("dofmap_edge_fine.mat",{"dofmap_edge_fine":dofmap_edge_fine})
savemat("dofmap_node_coarse.mat",{"dofmap_node_coarse":dofmap_node_coarse})
savemat("dofmap_node_fine.mat",{"dofmap_node_fine":dofmap_node_fine})

  




 
# v =Function(W)
# bmesh = BoundaryMesh(mesh2, 'exterior')
# coordinates_bmesh=bmesh.coordinates()


#from scipy.io import savemat
#savemat("mio.mat",Gmat.array())

#sigma1 = interpolate(sigma_raw1,RT)
#sigma2 = interpolate(sigma_raw2,RT)
#sigma_tensor= interpolate(sigma_tot,RTtensor)
#err_tensor=sigma_tensor - sigma
#sigmaxx = interpolate(sigmaxx_ex, feRT)
#sigmaxy = project(sigmaxy_ex, feRT)
#sigmayy = project(sigmayy_ex, feRT)

# Save solution in VTK format
#ufile_pvd = File("lsfem_elasticity_error_sigma.pvd")
#ufile_pvd << err_tensor




# Define Dirichlet boundary (x = 0 or x = 1)
#def boundary_L(x):
#    return x[0] < aL + DOLFIN_EPS
#def boundary_B(x):
#    return x[1] < bB + DOLFIN_EPS
#def boundary_R(x):
#    return x[0] > aR - DOLFIN_EPS
#def boundary_T(x):
#    return x[1] > bB - DOLFIN_EPS

#class BOUNDARY_TOP(SubDomain):
#	def inside(self,x,on_boundary):return ((x[1] > bB - DOLFIN_EPS) and on_boundary)
#class BOUNDARY_RIGHTBOTTOM(SubDomain):
#	def inside(self,x,on_boundary):return (((x[0] > aR - DOLFIN_EPS)or (x[1] < bB + DOLFIN_EPS)) and on_boundary)

# sigmaxx_ex=Expression("G *  pi * cos(pi*x[0])*sin(pi*x[1])+ L * pi * sin(pi*x[0])*cos(pi*x[1]) ", G=gamma,L=Lambda,degree=10)
# sigmaxy_ex=Expression("M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])",M=Mu,degree=4)
# sigmayy_ex=Expression("L * pi * cos(pi*x[0])*sin(pi*x[1])  + G  * pi * sin(pi*x[0])*cos(pi*x[1])",L=Lambda,G=gamma,degree=4)
# sigma_raw1=Expression(("G *  pi * cos(pi*x[0])*sin(pi*x[1])+ L * pi * sin(pi*x[0])*cos(pi*x[1]) ","M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])"), G=gamma,L=Lambda,M=Mu,degree=10)
# sigma_raw2=Expression(("M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])","L * pi * cos(pi*x[0])*sin(pi*x[1])  + G  * pi * sin(pi*x[0])*cos(pi*x[1])"), G=gamma,L=Lambda,M=Mu,degree=10)
# sigma_tot=Expression((("G *  pi * cos(pi*x[0])*sin(pi*x[1])+ L * pi * sin(pi*x[0])*cos(pi*x[1]) ","M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])"),("M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])","L * pi * cos(pi*x[0])*sin(pi*x[1])  + G  * pi * sin(pi*x[0])*cos(pi*x[1])")), G=gamma,L=Lambda,M=Mu,degree=10)
# h = Expression((("0.0","0.0"),("0.0","0.00")),degree=2)
#j = Expression(("1.0","0.00"),degree=2)
# j = Expression((("0.0","-0.005"),("0.0","0.0")),degree=2)
# kk = Expression(("0.0","-0.005"),degree=2)

#pc = mesh2.data().array("parent_cell", mesh2.topology().dim())
#print(pc)
# V1=FunctionSpace(mesh,"P",1)
# V2=FunctionSpace(mesh2,"P",1)
# mat=PETScDMCollection.create_transfer_matrix(V1,V2)
# print(mat.array())



#WC2WF.setValues([W_num_dofs-1],[W_num_dofs-1],1)
#WC2WF.assemble()
#print(WC2WF.getValues(range(15,W_num_dofs), range(15,W_num_dofs)))
#print(type(WC2WF))

#viewer=PETSc.Viewer().createBinary('test.dat', 'w')
#viewer(WC2WF)
#BB=WC2WF.convert("dense")
#BB.getDenseArray()
#print(BB)


#element = W.element()




#dofs_x = W.tabulate_dof_coordinates().reshape(W_num_dofs,gdim)


 
# for dofs, dofs_x in zip(dofs, dofs_x):
#      print(dofs)
#      print( ':') 
#      print(dofs_x)



# for c in cells(mesh):
#  print(W.dofmap().cell_dofs(c.index()))
  #    for v in vertices(c):
   #    print(v)
     #  Prova2[v]=Prova2[v]+1 





#marker = MeshFunction("bool", mesh, mesh.topology().dim(), True)




#print(dm2.entity_dofs(mesh2,0,np.array(nodes,dtype=np.uintp)))
#print(mesh_points)
#print(mesh_points2)
#print(dm.entity_dofs(mesh,0,np.array(nodes,dtype=np.uintp)))
#print(WC2WFnumpy)
#print()
#print("fine edges")

#for ee in range(0,30):
   #provami1=[Facet(mesh2,ee).midpoint().x(),Facet(mesh2,ee).midpoint().y()] 
   #provami2=dm2.entity_dofs(mesh2,1,np.array([ee],dtype=np.uintp))
   #print(ee)
   #print(provami1)
   #print(provami2)
    
#print()
#print("coarse edges")
#for ee in range(0,9):
   #provami1=[Facet(mesh,ee).midpoint().x(),Facet(mesh,ee).midpoint().y()] 
   #provami2=dm.entity_dofs(mesh,1,np.array([ee],dtype=np.uintp))
   #print(provami1)
   #print(provami2)
    
#print()
#print("fine nodes")
#for nn in range(0,15):
   #provami1=mesh_points2[nn]
   #provami2=dm2.entity_dofs(mesh2,0,np.array([nn],dtype=np.uintp))
   #print(provami1)
   #print(provami2)
#print()
#print("coarse nodes")
#for nn in range(0,6):
   #provami1=mesh_points[nn]
   #provami2=dm.entity_dofs(mesh,0,np.array([nn],dtype=np.uintp))
   #print(provami1)
   #print(provami2)
   
   
   
   ############## PETSC MATRIX
# WC2WF= PETSc.Mat().create()
# WC2WF.setSizes([W2_num_dofs, W_num_dofs])
# WC2WF.setType("aij")
# WC2WF.setUp()
   
# WC2WF= PETSc.Mat().create()
# WC2WF.setSizes([W2_num_dofs, W_num_dofs])
# WC2WF.setType("aij")
# WC2WF.setUp()
         #                WC2WF.setValues(tmp[0],[tmpC1[0], tmpC2[0]],[0.5, 0.5])
#                         WC2WF.setValues(tmp[1],[tmpC1[1], tmpC2[1]],[0.5, 0.5])

   #WC2WF.assemble()
   
   
   
   
   
   
   
   
   
   
   
   # boundaries = FacetFunction("size_t", mesh)
# boundaries.set_all(0)
# BOUNDARY_TOP().mark(boundaries, 1)
# BOUNDARY_RIGHTBOTTOM().mark(boundaries, 2)

# ds = Measure('ds', domain = mesh, subdomain_data = boundaries)

#J =C_constitutive * inner(A(sigma)-epsilon(u),A(sigma)-epsilon(u))*dx + C_equilibrium * inner(div(sigma)+f,div(sigma)+f)*dx+ C_boundary* inner(dot(sigma,FacetNormal(mesh))-mm,dot(sigma,FacetNormal(mesh))-mm)*ds(1)+C_boundary* inner(dot(sigma,FacetNormal(mesh))-nn,dot(sigma,FacetNormal(mesh))-nn)*ds(2)

#bc1 = DirichletBC(W.sub(1), g, boundary_L)
#bc2 = DirichletBC(W.sub(1), g, boundary_B)
#bc4 = DirichletBC(W.sub(1), g, boundary_T)
#bc3 = DirichletBC(W.sub(1), g, boundary_R)


#bc1 = DirichletBC(W.sub(1), g, boundary_L)
#bc2 = DirichletBC(W.sub(0), h, boundary_B)
#bc4 = DirichletBC(W.sub(0).sub(1), kk, boundary_T)
#bc3 = DirichletBC(W.sub(0), h, boundary_R)
#bcs=[bc1,bc2,bc4,bc3]



# w.vector()[:] = 1 #0
# w.vector()[1]=0
# w.vector()[2]=0
# w.vector()[12]=0
# w.vector()[14]=0
# w.vector()[13]=0
# w.vector()[15]=0
#constant_vector = Constant(("1.0","1.0","1.0","1.0","1.0","1.0"))
#w = interpolate(constant_vector, W)

