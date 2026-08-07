import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.ErsatzAnwendung import *
kernel = Kernel()   #defining kernel

data = [[2, 4], [1, 3], [0, 0], [0, 0]]
#REF http://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm

A = Matrix(4, 2)
for i in range(0, 4):
    for j in range(0, 2):
        A[i, j] = data[i][j]
print("A:" + str(A))

util = POD_Utils()
Out = util.SVD(A)
print("U:" + str(Out[0]))
print("S:" + str(Out[1]))
print("VT:" + str(Out[2]))

