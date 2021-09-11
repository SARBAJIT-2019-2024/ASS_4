#         #   #     #         #            #     # # # #        # # # #      #     #
# #     # #    #   #          #            #     #     #        #     #       #   #
#  #   #  #     # #           #            #     #     #        #     #        # #
#   # #   #      #            #            #     # # # #        # # # #         #
#    #    #      #            #            #     #     #        #               #
#         #      #            #            #     #     #        #               #  
#         #      #   #######  # # # # #    #     # # # #   #    #               # 

#####################################################
#Calling a matrix
''''''
#list_C=[]
#with open("txt file") as matC:
    #for k in matC:
        #list_C.append(list(map(float, k.split()))

       
                

## Custom library for importing functions
##################################################################################################
#Partial Pivot
def partial_pivoting(a,b,n):
    n=len(a)
    no_of_swaps=0
    for i in range(n-1):
        if abs(a[i][i]) == 0:
            for j in range(i+1,n):
                if abs(a[j][i]) > abs(a[i][i]):
                    a[j], a[i] = a[i], a[j]  # interchange ith and jth rows of matrix 'A'
                    b[j], b[i] = b[i], b[j]  # interchange ith and jth elements of vector 'b'
            no_of_swaps+=1
    return no_of_swaps        

###################################################################################################
#defining Gauss Jordan
def Gauss_jordan(list_C):
    a,b=Matrixmaker(list_C)
    n=len(a)
    
    partial_pivoting(a,b,n) #calling the partial pivot
    
    for i in range(n):
        pivot_element = a[i][i]
        
        
        if pivot_element==0:      #condition if pivot element turns out be zero in any step
            return (None,None)
        if type(b[i]) is list: 
            for l in range(n):
                if b[i][l] != 0:
                    b[i][l] = b[i][l]/pivot_element
            for r in range(i, n):
                a[i][r] = a[i][r]/pivot_element
            for k in range(n):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    for j in range(i,n):
                        a[k][j] = a[k][j] - balance_factor*a[i][j]  
                    for l in range(n):
                        if b[i][l] != 0:
                            b[k][l] = b[k][l] - balance_factor*b[i][l]

                        
        else: 
            b[i] = b[i]/pivot_element                  #Condition for pivot element
            for k in range(i, n):
                 a[i][k] = a[i][k]/pivot_element       #dividing the elements by pivot element    

            for k in range(n):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    b[k] = b[k] - balance_factor*b[i]
                    for j in range(i, n):
                        a[k][j] = a[k][j] - balance_factor*a[i][j] 
    return(a,b)
###############################################################################################
# Defining matrix maker which makes a matrix C in to A and b
    
def Matrixmaker(list_C):
    list_A=[[0 for x in range(len(list_C))] for y in range(len(list_C))]
    for i in range(len(list_C)):
            for j in range(len(list_C)):
                list_A[i][j]=list_C[i][j] 
    if len(list_C[0])==len(list_C)+1:          #Generating list of C  
        list_B=[0 for x in range(len(list_C))] #generating a list of B where we store the result
        for i in range(len(list_C)):
            list_B.append(0)
        
        for i in range(len(list_C)):
            list_B[i]=list_C[i][len(list_C)]
    else:
        list_B=[[0 for x in range(len(list_C))] for y in range(len(list_C))] 
        for i in range(len(list_C)):

            list_B[i][i]=1       
    return(list_A,list_B)

#####################################################################################################
 #Defining matrix multiplication

def matrix_mul(a,b):
    AB_=[[0 for x in range(len(a))] for y in range(len(a))] 
        
    for i in range(len(a)):
        for j in range(len(b[2])):#expanding along the coloumn 2 of B matrix
            for k in range(len(b)):
                AB_[i][j] += a[i][k]*b[k][j] 
    return AB_
#####################################################################################################    

# Defining function for determinant calculation

def determinant_calc(a):           
    if len(a) != len(a[1]):         
        print("The determinant is not defined, Provide  Square Matrix please !")#checking for entry of Square matrix.
    else:     
        no_of_swaps = 0
        for i in range(len(a)): 
            if abs(a[i][i])== 0:
                for j in range(i+1 , len(a)): 
                    if  abs(a[i][i]) < abs(a[j][i]): 
                        for k in range(i , len(a)):
                            a[i][k], a[j][k] = a[j][k], a[i][k]         #swapping the rows of the matrix.
                no_of_swaps += 1  
        for i in range(len(a)):
            pivot_element = a[i][i]
            if pivot_element==0:
                return (0)
            for k in range(len(a)):
                if k != i and a[k][i] != 0:
                    balance_factor = a[k][i]
                    for j in range(i,len(a)):
                        a[k][j] = a[k][j] - balance_factor*(a[i][j]/pivot_element) #making the elements zero by balance factor
        det = 1
        for i in range(len(a)):
            det*=a[i][i]
        det*=(-1)**(no_of_swaps)    
        print("The determinant for the given matrix is:")        
            
        return(det)     
                
#################################################################################################################################
#################################################################################################################################







#CROUT METHOD

def crout_method(a):
    n = len(a)
    for i in range(n):
        for j in range(i,n):
            sum_lower = 0
            for k in range(i):
                sum_lower += (a[i][k]*a[k][j])
            a[i][j] = a[i][j] - sum_lower
        for j in range(i,n):
            if (i == j):
                pass
            else:
                sum_upper = 0
                for k in range(i):
                    sum_upper += (a[j][k]*a[k][i])
                a[j][i] = (a[j][i]- sum_upper)/a[i][i]
    Lower_tri, Upper_tri = [[0 for x in range(n)] for y in range(n)], [[0 for x in range(n)] for y in range(n)] #Breaking the matrix into two parts upper and lower
    partial_pivoting(Upper_tri,Lower_tri,n)
    for i in range(n):
        for j in range(i+1):
            if (i == j):     # conditioning for crout U[i][i]=1
                Upper_tri[j][j] = 1
                Lower_tri[j][j] = a[j][j]
            else:
                Upper_tri[j][i], Lower_tri[i][j] = a[j][i], a[i][j]
    return  Upper_tri,Lower_tri
            
#####################################################################################################################################
#####################################################################################################################################


def combiner(Upper_tri,Lower_tri):  # combining upper and lower matrix in one matrix 
    n=len(Upper_tri)
    D=[[0 for x in range(n)] for y in range(n)]
    for i in range(len(Upper_tri)):
        for j in range(len(Upper_tri)):
            if i==j:
                D[j][j]=Lower_tri[i][j]   #replacing diagonal elements by Lower triangular's diagonal elements
            else:
                D[i][j]=Lower_tri[i][j]+Upper_tri[i][j]
    return D            


   
####################################################################################################################################
####################################################################################################################################

#FORWARD BACKWARD SUB FOR DO LITTLE

def forward_backward_sub_do_little(Upper_tri, Lower_tri, b):
    n=len(Upper_tri)
    y = [0 for i in range(n)]

    for i in range(n):
        total = 0
        for j in range(i):
            total += Lower_tri[i][j] * y[j]
        y[i] = b[i] - total

    x = [0 for i in range(n)]

    for i in reversed(range(n)):
        total = 0
        for j in range(n):
            total += Upper_tri[i][j] * x[j]
        
        x[i] = (y[i] - total)/Upper_tri[i][i]

    return x
########################################################################################################################
########################################################################################################################

#FORWARD BACKWARD SUB FOR CROUT

def forward_backward_sub_crout(Upper_tri, Lower_tri, b):
    n=len(Upper_tri)
    partial_pivoting(Upper_tri,Lower_tri,n)
    D= combiner(Upper_tri,Lower_tri)
    
    y = [0 for i in range(n)]

    for i in range(n):
        total = 0
        for j in range(i):
            total += D[i][j] * y[j]
        y[i] = b[i] - total

    x = [0 for i in range(n)]

    for i in reversed(range(n)):
        total = 0
        for j in range(n):
            total += D[i][j] * x[j]
        
        x[i] = (y[i] - total)/D[i][i]

    return x


  
#####################################################################################################################
#####################################################################################################################

#CROUT METHOD FOR solution
def linear_solver_crout(list_C):
    A,b=Matrixmaker(list_C)
    n= len(A)
    
    partial_pivoting(A, b, n)
    for i in range(len(A)):
        if A[i][i]==0:
            print("No solutions exist for this system !")
            return None
        
    Upper_tri, Lower_tri = crout_method(A)
    x = forward_backward_sub_crout(Upper_tri,Lower_tri, b)
    print("The solutions of the system of linear equations by Crout's method is")
    return x

################################################################################################################
################################################################################################################

# UPPER LOWER BY DO-LITTLE
def Do_little(a):
    
    n= len(a)
    Lower_tri = [[0 for x in range(n)]
             for y in range(n)]
    Upper_tri = [[0 for x in range(n)]
             for y in range(n)]
 
    for i in range(n):
        Lower_tri[i][i] = 1

    for j in range(n):
        for i in range(n):
            total = 0
            for k in range(i):
                total +=Lower_tri[i][k] * Upper_tri[k][j]

            if i == j:
                Upper_tri[i][j] = a[i][j] - total

            elif i > j:
                Lower_tri[i][j] = (a[i][j] - total)/Upper_tri[j][j]

            else :
                Upper_tri[i][j] = a[i][j] - total

    return Upper_tri, Lower_tri
################################################################################################################
################################################################################################################
#Dolittle METHOD FOR solution
def linear_solver_do_little(list_C):
    A,b=Matrixmaker(list_C)
    n= len(A)
    
    partial_pivoting(A, b, n)
    for i in range(len(A)):
        if A[i][i]==0:
            print("No solutions exist!")
            return None
    Upper_tri, Lower_tri = Do_little(A)
    x = forward_backward_sub_do_little(Upper_tri,Lower_tri, b)
    print("The solutions of the system of linear equations by Dolittle's method is")
    return x

##################################################################################################################
##################################################################################################################

#printing upper and lower triangular matrix

def print_matrix(Lower_tri,Upper_tri):
    
    for i in range(len(Lower_tri)):
        if Lower_tri[i][i]==1:
            for j in range(len(Lower_tri)):
                if i>=j:
                    print(Upper_tri[i][j],end=" ")
                else:
                    print(Lower_tri[i][j],end=" ")  
            print()          
        
        else :
            for j in range(len(Lower_tri)):
                if i>j:
                    print(Upper_tri[i][j],end=" ")
                else:
                    print(Lower_tri[i][j],end=" ")  
            print() 
    return  None

#####################################################################################################################
#####################################################################################################################

def Nullmatrix(m,n): # creating a null matrix
	p= [[0 for i in range(n)] for j in range(m)]
	return(p)
#Creating square matrix
def identity_mat(m):
	p=Nullmatrix(m,m)
	for i in range(m):
		p[i][i] = 1
	return(p)

###################################################################################################################
###################################################################################################################
def pivoting_matrix(a): #PA=LU
    
    m = len(a)

    # Creating an identity matrix                                                                                                                                                                                          
    identity_matrix = [[float(i ==j) for i in range(m)] for j in range(m)]

    #The largest element in each row will be placed into the diagonal                                                                                                                                                                                             
    for j in range(m):
        row = max(range(j, m), key=lambda i: abs(a[i][j]))
        if j != row:
            # Swapping the rows                                                                                                                                                                                                                            
            identity_matrix[j], identity_matrix[row] = identity_matrix[row], identity_matrix[j]

    return identity_matrix
################################################################################################################################
################################################################################################################################

#inverse using LU inverse
def LU_inverse(a):
    
    
    
    I=identity_mat(4)
    P = pivoting_matrix(a) #calling the pivoting matrix
    PA = matrix_mul(P, a) #PA=LU

    P_1=Nullmatrix(4,4) # as P_inverse= P_transpose 
    for i in range(4):
	    for j in range(4):
		    P_1[j][i]=round(P[i][j],3)


    Upper_tri,Lower_tri=crout_method(PA)
    
    det=1            # condition for inverse existing
    for i in range(len(Lower_tri)):
        det*=Lower_tri[i][i]
    if det==0:
        print("Inverse doesn't exist!")
        return None     
# stroing matrix
    x_1=Nullmatrix(4,4)
    x_2=Nullmatrix(4,4)  


    for i in range(4):
	    x_1[i]=forward_backward_sub_crout(Upper_tri,Lower_tri,I[i])

    Px_i = matrix_mul(P_1,x_1)	


    for i in range(4):
	    for j in range(4):
		    x_2[j][i]=round(Px_i[i][j],5)

    return x_2



################################################################################################################
################################################################################################################    

#Forward backward method for cholesky's substitution 


def forward_backward_sub_cholesky(Lower_tri, Lower_tri_Transpose, b):
    n=len(Lower_tri)
    
    
    
    y = [0 for i in range(n)]

    for i in range(n):
        total = 0
        for j in range(i):
            total += Lower_tri[i][j] * y[j]
           
        y[i] = (b[i] - total)/Lower_tri[i][i]

    x = [0 for i in range(n)]

    for i in reversed(range(n)):
        total = 0
        for j in range(n):
            total += Lower_tri_Transpose[i][j] * x[j]
        
        x[i] = (y[i] - total)/Lower_tri_Transpose[i][i]

    return x


#################################################################################################################
#################################################################################################################

#CHOLESKY's METHOD

def cholesky_decomp (C):
    a,b=Matrixmaker(C)
    n=len(a)
    Lower_tri = [[0 for x in range(n)]
             for y in range(n)]
    Lower_tri_Transpose = [[0 for x in range(n)]
             for y in range(n)]         
    

    for i in range(n):
        for j in range(i + 1):
            total = 0
 
            if j == i:
                for k in range(j):
                    total +=Lower_tri[i][k] **2
            
                Lower_tri[j][j] =(a[i][i] - total)**(0.5)
                
            else:
                for k in range(j):
                    total += (Lower_tri[i][k] *Lower_tri[j][k])
                if(Lower_tri[j][j] > 0):
                    Lower_tri[i][j] = ((a[i][j] - total) /
                                               Lower_tri[j][j])
    for i in range(len(Lower_tri)):
        for j in range(len(Lower_tri)):           
                Lower_tri_Transpose[i][j]=Lower_tri[j][i]
    
    return Lower_tri,Lower_tri_Transpose

def Cholesky_Solver(A):
    a,b=Matrixmaker(A)
    Lower_tri,Lower_tri_Transpose=cholesky_decomp(a)
    x = forward_backward_sub_cholesky(Lower_tri,Lower_tri_Transpose, b)
    print("The solutions of the system of linear equations by Cholesky's method is")
    return x

##################################################################################################################
##################################################################################################################


