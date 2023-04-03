import numpy as np
np.set_printoptions(precision=7, suppress=True, linewidth=100)


def og_ode(t:float, y:float):
  return t - (y ** 2)

def ode_first_deriv(t:float, y:float):
  return og_ode(t,y) - 2 * y





def do_shit(t, y, h):
  yi = og_ode(t,y)
  hy_ti = h * ode_first_deriv(t,y)
  return yi + hy_ti
  

def eulers():
  f0 = 1
  num_iter = 10
  start, end = (0,2)

  h = (start - end)/num_iter

  for cur_iter in range(0,10):
    t =  cur_iter
    y = f0
    h = h
    inner_math = do_shit(t,y,h)
    next_f0 = inner_math
    
    
    f0 = next_f0
    cur_iter 
    
    
  return next_f0

def runge_kutta():
  f0 = 1
  num_iter = 10
  start, end = (0,2)
  ti = 0
  yi = 1

  h = 0.2
  for i in range(0,10):
    k1 = ode_first_deriv(ti, yi)
    k2 = ode_first_deriv(ti + h/2, yi + (h/2) * k1)
    k3 = ode_first_deriv(ti + h/2, yi + (h/2) * k2)
    k4 = ode_first_deriv(ti + h, yi + h * k3)
    yi_next = yi + (h/6) * (k1 + (2 * k2) + (2 * k3) + k4)
    
    yi = yi_next
    i = i + 1
    ti = ti + h
  return yi_next



      

  



def gje_row_ops(A,b):
  n = len(b)
  Ab = np.concatenate((A, b.reshape(n,1)), axis=1)
 
  Ab[0] = Ab[0]*0.5
  
  Ab[2] = -1 * Ab[2]
  
  Ab[1] = Ab[1] - Ab[0]
 
  Ab[2] = Ab[2] - Ab[0]
  
  Ab[1] = (2*Ab[1])/7
  
  
  Ab[2] = (-2*Ab[2])/9
  
  
  Ab[2] = Ab[2] - Ab[1]
  
  Ab[2] = (7*Ab[2])/6
  
  
  Ab[1] = Ab[1] - (Ab[2]/7)
  
  Ab[0] = Ab[0] - (Ab[2]*0.5)
  
  Ab[0] = Ab[0] + (Ab[1]*0.5)
  
  result = [Ab[0,3],Ab[1,3],Ab[2,3]]
  return result

def LU_row_ops(A):
  L = np.array([[1.0,0.0,0.0,0.0],
                [0.0,1.0,0.0,0.0],
                [0.0,0.0,1.0,0.0],
                [0.0,0.0,0.0,1.0]])
  U = A

  U[1] = U[1] - (2*U[0])
  L[1,0] = 2.0

  U[2] = U[2] - (3*U[0])
  L[2,0] = 3.0

  U[3] = U[3] - (-U[0])
  L[3,0] = -1.0

  U[2] = U[2] - (4*U[1])
  L[2,1] = 4.0

  U[3] = U[3] - (-3*U[1])
  L[3,1] = -3.0

  U[3] = U[3] - (0*U[2])
  L[3,2] = 0.0

  L = np.array([[1.0,0.0,0.0,0.0],
               [2.0,1.0,0.0,0.0],
               [3.0,-4.0,1.0,0.0],
               [-1.0,-3.0,0.0, 1.0]])
  detA = U[0,0] * U[1,1] * U[2,2] * U[3,3]
  print(detA)
  print("\n")
  print(L)
  print("\n")
  print(U)

def is_diag_dom(A):
  sum1 = A[0,0] + A[0,1] + A[0,2] + A[0,3] + A[0,4]
  sum2 = A[1,0] + A[1,1] + A[1,2] + A[1,3] + A[1,4]
  sum3 = A[2,0] + A[2,1] + A[2,2] + A[2,3] + A[2,4]
  sum4 = A[3,0] + A[3,1] + A[3,2] + A[3,3] + A[3,4]
  sum5 = A[4,0] + A[4,1] + A[4,2] + A[4,3] + A[4,4]
  if((A[0,0] >= sum1) & (A[1,1] >= sum2) & (A[2,2] >= sum3) & (A[3,3] >= sum4) & (A[4,4] >= sum5)):
    dd = True
  else:
    dd = False
  return dd

def is_pos_def(A):
  row_len = len(A[0])
  
  col_len = len(A[:,0])

  if(row_len == col_len):
    pd = True
  else:
    pd = False
  return pd
  
if __name__ == "__main__":
  print(eulers())
  print("\n")
  print(runge_kutta())
  print("\n")
  A1 = np.array([[2.0,-1.0,1.0],
                [1.0,3.0,1.0],
                [-1.0,5.0,4.0]])
  b = (np.array([6.0,0.0,-3.0]))
  x = gje_row_ops(A1,b)
  print(x)
  print("\n")
  A2 = (np.array([[1.0,1.0,0.0,3.0],
                  [2.0,1.0,-1.0,1.0],
                  [3.0,-1.0,-1.0,2.0],
                  [-1.0,2.0,3.0,-1.0]]))
  
  LU_row_ops(A2)
  print("\n")
  A3 = np.array([[9.0,0.0,5.0,2.0,1.0],
                [3.0,9.0,1.0,2.0,1.0],
                [0.0,1.0,7.0,2.0,3.0],
                [4.0,2.0,3.0,12.0,2.0],
                [3.0,2.0,4.0,0.0,0.8]])
  print(is_diag_dom(A3))
  print("\n")
  A4 = np.array([[2,2,1],
                [2,3,0],
                [1,0,2]])
  print(is_pos_def(A4))

