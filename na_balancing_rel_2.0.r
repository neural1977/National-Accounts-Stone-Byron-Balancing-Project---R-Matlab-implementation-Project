###################################################################################
#                                                                                 #
#         National Accounts Balancing Procedure                                   #            
#                                                                                 #
#         Implementation: R Language                                              #
#         Date of creation: 06/06/2013                                            #
#         Date of last change: 31/07/2013                                         #
#         Release Version: 2.0                                                    #                          
#                                                                                 #
#         Author: Francesco Pugliese                                              #
#         Email: pugliese05@gmail.com                                             #                          
#                                                                                 #
###################################################################################

##########################################
#     General Direct Stone Balancing     #
##########################################

########################
#   GLOBAL SETTINGS    #
########################

####   SELECTION OF METHOD SECTION  (DIRECT, ITERATIVE, ETC..) ####
direct=2                                                                              # 0=direct method - Balancing Method by Stone's Estimator,1=iterative method - Balancing Method by Conjugate Gradient Method, 2=iterative method CG with pre-conditioning matrix W
methods_names=c("Direct Method - Stone's Estimator","Iterative Method - Basic Conjugate Gradient","Iterative Method - Conjugate Gradient by Pre-conditioning Matrix \"W\"")
number_of_implemented_methods=3                                                       # Number of all implemented methods
test_all=0                                                                            # 1= test all methods, 0=test only the specified one

####   PATH SETTINGS    ####
#experiment_name="costruzioni"                                                        # name of the balancing experiment
#experiment_name="simulation_test"                                                    # name of the balancing experiment
experiment_name="global_balancing_2"                                                  # name of the balancing experiment
data_input="C:/Users/frpuglie/Desktop/Stone Balancing/Bilanciamento in R/input"	      # Data input, personalize this folder according the needs
data_output="C:/Users/frpuglie/Desktop/Stone Balancing/Bilanciamento in R/results"	  # Data output, personalize this folder according the needs


####   SIMULATION SETTINGS   ####
simulation_test=0                                                                     # Enable a loop of testing with random matrix
type_of_simulation=1                                                                  # Type of simulation: 1-tests on matrix size 2-tests on matrix density
number_of_tests=1                                                                     #=1, only to test on one random matrix with specified dimensions
max_matrix_row_size=22
max_matrix_col_size=22
max_matrix_density=100                                                   
max_random_matrix_element_values=10000                                                # Max of generated random matrix values  
max_random_marginals_perturbation=200                                                 # Max of perturbation in random marginals  

####   FILES LOADING SETTINGS    ####
save_to_xls=0                                                                         # Save output matrices to xls 
load_from_xls=0                                                                       # Loads matrices from xls and deactivate simulations, etc.. 1=load from file, 0=run the built-in example   
v_matrix=1                                                                            # 0=c matrix given by input, diag(v)=c*p, 1=v matrix, already given by input  
input_matrix.filename="matrix.xls"                                                    # Name of the input matrix file
input_row_marginals.filename="offert.xls"                                             # Name of the input row marginals file
input_col_marginals.filename="domand.xls"                                             # Name of the input col marginals file
v_matrix.filename="vmatrix.xls"                                                       # Name of the input V matrix file
v_row_marginals.filename="voffert.xls"                                                # Name of the input V row marginals file
v_col_marginals.filename="vdomand.xls"                                                # Name of the input V col marginals file
c_matrix.filename=""                                                                  # Name of the input C matrix file
c_row_marginals.filename=""                                                           # Name of the input C row marginals file
c_col_marginals.filename=""                                                           # Name of the input C col marginals file

####   PARSING CONFIGURATION FILE   ####
load_from_deck=1                                                                      # Load from configuration file
config_file_name="vincoli.deck"
config_language_re="([-,+]|MM|SR|SC|VR|VC)"                                           # Regular expression of configuration file's language
sub_domains=0                                                                         # Generates balancing sub-domains

####  INNER TRIAL EXAMPLE SETTINGS   ####
matrix_row_size=3                                                                     # matrix row size
matrix_col_size=3                                                                     # matrix column size


####   FILES LOADING SETTINGS    ####
plot_all=1                                                                            # Plot all the related charts  
save_plot_to_file=0                                                                   # Save the charts on a file: 1=jpeg, 2=pdf
plot_type=1                                                                           # type of plot: 1=line plot, 0=bar plot

####   MISCELLANEOUS    ####
sep_windows=1		                                     									                # 0=no sep windoes, 1=sep windows
lim_perc=15 											                                                    # Percentage of plot scaling
compress=1                                                                            # 1=compress sparse matrices to enhance times and memory usage     
exit_point=1e-10                                                                      # Exit error in conjugate gradient method

if (compress)
    library(Matrix) 
        
############
#  LOADER  #
############

#ptm <- proc.time()                                                                    # Time acquisition 

if (load_from_xls) {
  library(xlsx)
  simulation_test=0                                                                                                                             # deactivate simulation_test
  xls_in_matrix=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",input_matrix.filename,sep=""),1, header=F))                        # read input matrix from excel file, first sheet
  xls_in_row_marginals=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",input_row_marginals.filename,sep=""),1, header=F))          # read input col marginals from excel file, first sheet
  xls_in_col_marginals=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",input_col_marginals.filename,sep=""),1, header=F))          # read input row marginals from excel file, first sheet
  
  if (v_matrix) {
    xls_v_matrix=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",v_matrix.filename,sep=""),1, header=F))                             # read v matrix from excel file, first sheet
    xls_v_row_marginals=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",v_row_marginals.filename,sep=""),1, header=F))               # read v col marginals from excel file, first sheet
    xls_v_col_marginals=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",v_col_marginals.filename,sep=""),1, header=F))               # read v row marginals from excel file, first sheet
  } else {
    xls_c_matrix=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",c_matrix.filename,sep=""),1, header=F))                             # read c matrix from excel file, first sheet
    xls_c_row_marginals=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",c_row_marginals.filename,sep=""),1, header=F))               # read c col marginals from excel file, first sheet
    xls_c_col_marginals=as.matrix(read.xlsx(paste(data_input,"/",experiment_name,"/",c_col_marginals.filename,sep=""),1, header=F))               # read c row marginals from excel file, first sheet
  }

  if (experiment_name=="costruzioni") {
    xls_in_matrix=xls_in_matrix[,-c(3,4,10,12,14,15,18)]  
    xls_in_col_marginals=xls_in_col_marginals[,-c(3,4,10,12,14,15,18)]
    xls_v_matrix=xls_v_matrix[,-c(3,4,10,12,14,15,18)]  
    xls_v_col_marginals=xls_v_col_marginals[,-c(3,4,10,12,14,15,18)]
  }

  #Loaded matrix dimension, define new ones
  matrix_row_size=dim(xls_in_matrix)[1]                                       
  matrix_col_size=dim(xls_in_matrix)[2]                   
}

##############
#  SETTINGS  #
##############

if (save_plot_to_file)                                                                #note: set sep_window to 0 if save_plot_to_file=1
  sep_windows=0

# Parsing of configuration file
if (load_from_deck) {
  simulation_test=0
  load_from_xls=0
  
  all_data = readLines(paste(data_input,"/",experiment_name,"/", config_file_name, sep=""))       # load all the lines
  data_size=dim(as.matrix(all_data))[1]
  
  if (sub_domains) {                                                                              # parse and generates balancing sub_domains to simplify the whole problem
    whole_set=1:data_size
  
    # Section for the Reduction of the balancing problem to k consistent balancing sub-domains
    if (compress==0) {
      consistent_blocks=array(0,dim=c(data_size, data_size))                                      # size the consistent of blocks
      tmp_consistent_block=array(0,dim=data_size*data_size)                                       # size the temporary general consistent of blocks
      tmp_unique=array(0,dim=data_size)                                                           # size the temporary consistent of blocks
    } else {
      consistent_blocks=Matrix(nrow = data_size, ncol = data_size, 0, sparse = TRUE)              # size the consistent of blocks in optimized form 
      tmp_consistent_block=Matrix(nrow = data_size*data_size, 0, sparse=TRUE)                     # size the temporary general consistent of blocks in optimized form
      tmp_unique=Matrix(nrow = data_size, 0, sparse=TRUE)                                         # size the temporary consistent of blocks in optimized form
    }

    # Main loop through the sub-blocks
    block_number=1
    #for (block_number in 1:2) {                                                                       # loop on block number  
      current_index=0
      lll=1                                                                                           # initialize whole_set index
      l=whole_set[lll]                                                                                  # first line
      ll=whole_set[lll]                                                                                 # first line
    
      # Loop through the sub-block's lines 
      while(lll<data_size+1) {  
        prev_index=0
        current_line=strsplit(all_data[l],"\\s+")[[1]]                                                # parse the single line in order to find the objects
        token_set=grep(config_language_re, current_line, invert=TRUE)                                 # determines all the tokens to search in the file
        ntokens=dim(as.array(token_set))                                                              # number of tokens
  
        #build a balancing consistent sub-block
        for (tt in 1:ntokens) {
          token=current_line[token_set[tt]]                                                           # tokenize each line
          rep_data_entries = grep(token, all_data)                                                    # search other occurencies of the read token
          tmp_unique[(prev_index+1):(dim(as.matrix(rep_data_entries))[1]+prev_index)]=rep_data_entries           # build a block of indexes con consistent balancing elements
          prev_index=prev_index+dim(as.matrix(rep_data_entries))[1]
        }
        #removes duplicates
        tmp=unique(as.array(tmp_unique))                                                              # removes doubles 
        tmp=sort(tmp[-grep("([^0])", tmp, invert=TRUE)])                                              # removes all the 0 by the last grep        

        tmp_consistent_block[((l-1)*current_index+1):(((l-1)*current_index)+dim(as.array(tmp)))]=tmp  # place the new temporary array-block in the general block
        tmp=unique(as.array(tmp_consistent_block))                                                    # removes doubles 
        tmp=sort(tmp[-grep("([^0])", tmp, invert=TRUE)])                                              # removes all the 0 by the last grep        
    
        if (compress==0)
          tmp_consistent_block=array(0,dim=data_size*data_size)                                       # size the temporary general consistent of blocks
        else
          tmp_consistent_block=Matrix(nrow = data_size*data_size, 0, sparse=TRUE)                     # size the temporary general consistent of blocks in optimized form
    
        current_index=dim(as.array(tmp))+1                                                            # save the size of tmp for the next loop
        tmp_consistent_block[1:(current_index-1)]=tmp                                                 # place the new temporary array-block in the general block
  
        lll=lll+1                                                                                     # increases the whole_set index
        ll=whole_set[lll]                                                                              # browse inside whole_set
        l=tmp_consistent_block[ll]                                                                    # find the next index to browse
        if (l==0)
          break
      }
    
      tmp_consistent_block[1:(current_index-1)]=tmp
      consistent_blocks[1:(current_index-1),block_number]=tmp_consistent_block[1:(current_index-1)]    # load the found sub-block into the structure 
  
      whole_set=setdiff(whole_set,consistent_blocks[,block_number])                                    # new set under examination
      data_size=dim(as.matrix(whole_set))[1]                                                           # size of the new set
    #}
  } else {                                                                                             # simply parse the input configuration file
  
  
  }
}
                               
if (simulation_test) {
  if (type_of_simulation==1) {                                                        # tests on matrix size                                                           
    matrix_row_step=round(max_matrix_row_size/number_of_tests)                        # this rounds the size
    matrix_col_step=round(max_matrix_col_size/number_of_tests)
    matrix_density=max_matrix_density                                                   
  } else if (type_of_simulation==2) {                                                 # tests on matrix density
    matrix_row_step=max_matrix_row_size                                               
    matrix_col_step=max_matrix_col_size
    matrix_density_step=round(max_matrix_density/number_of_tests)                     # Matrix Density Step in percentage
  }
} else {
  number_of_tests=1                                                                   # in the single case, matrix size is the default one
  matrix_row_step=matrix_row_size
  matrix_col_step=matrix_col_size
}

#selection of the loop on methods
if (test_all) {
  number_of_methods=number_of_implemented_methods
} else {
  number_of_methods=1
} 

#Additional arrays
times_array=array(0,dim=c(5,number_of_tests,number_of_methods))                       # array of times 
errors_array=array(0,dim=c(4,number_of_tests,number_of_methods))                      # array of errors
if (!plot_type) {
  times_tmp_array=array(0,dim=c(number_of_methods, number_of_tests))                  # array of times temporary
  iters_tmp_array=array(0,dim=c(number_of_methods, number_of_tests))                  # array of iterations temporary
}

#Build colors array for plots
line_colors=array(0,dim=number_of_methods)
sym_table=array(0,dim=number_of_methods)
line_colors[1]="black"                                                                # colors of lines and legends
sym_table=21:(21+number_of_methods-1)                                                 # symbols of points in the legend  
# eventually adds other methods' data
if (number_of_methods>1)
  for (m in 2:number_of_methods)
    line_colors[m]=grey(0:(2*number_of_methods)/(2*number_of_methods))[2*m]           # generates grey tones

for (m in 1:number_of_methods) {
  if (test_all)
    direct=m-1                                                                        # loop through all the methods  
  
  for (g in 1:number_of_tests) {

    if (direct==0)
      print("DIRECT METHOD INITIALIZATION")
    else if (direct==1)
      print("CONJUGATE GRADIENT METHOD INITIALIZATION")
    else
      print("CONJUGATE GRADIENT METHOD WITH W INITIALIZATION")
  
    ptm <- proc.time()                                                                # Time acquisition 
 
    if (simulation_test) {
      if (type_of_simulation==1) {                                                    # tests on matrix size                                                           
        in_matrix_row_size=g*matrix_row_step
        in_matrix_col_size=g*matrix_col_step
      } else if (type_of_simulation==2) {                                             # tests on matrix density
        in_matrix_row_size=matrix_row_step
        in_matrix_col_size=matrix_col_step
        matrix_density=g*matrix_density_step                                          # increase the matrix density
      }
    } else {
      in_matrix_row_size=matrix_row_step
      in_matrix_col_size=matrix_col_step
    }
  
    ############
    #   MAIN   #
    ############

    # Declarations #
    if (compress==0) {
      in_matrix=array(0,dim=c(in_matrix_row_size, in_matrix_col_size))                            # allocates the input matrix for the balancing procedure
      row_marginals=t(array(0,dim=in_matrix_col_size))                                            # allocates the input column vector of row marginals
      col_marginals=array(0,dim=in_matrix_row_size)                                               # allocates the input row vector of col marginals

      c_matrix=array(0,dim=c(in_matrix_row_size, in_matrix_col_size))                             # allocates the reliability matrix for the balancing procedure
      c_row_marginals=t(array(0,dim=in_matrix_col_size))                                          # allocates the reliability column vector of row marginals
      c_col_marginals=array(0,dim=in_matrix_row_size)                                             # allocates the reliability row vector of col marginals
    } else { 
      in_matrix=Matrix(nrow = in_matrix_row_size, ncol = in_matrix_col_size, 0, sparse = TRUE)    # allocates the input matrix for the balancing procedure in optimized form
      row_marginals=t(Matrix(ncol = in_matrix_col_size, 0, sparse = TRUE))                        # allocates the input column vector of row marginals in optimized form
      col_marginals=Matrix(nrow = in_matrix_row_size, 0, sparse = TRUE)                           # allocates the input row vector of col marginals in optimized form

      c_matrix=Matrix(nrow = in_matrix_row_size, ncol = in_matrix_col_size, 0, sparse = TRUE)     # allocates the reliability matrix for the balancing procedure in optimized form
      c_row_marginals=t(Matrix(ncol = in_matrix_col_size, 0, sparse = TRUE))                      # allocates the reliability column vector of row marginals in optimized form
      c_col_marginals=Matrix(nrow = in_matrix_row_size, 0, sparse = TRUE)                         # allocates the reliability row vector of col marginals in optimized form
    }

    n=in_matrix_row_size*in_matrix_col_size+in_matrix_row_size+in_matrix_col_size        # n=p vector size
    k=in_matrix_row_size+in_matrix_col_size                                              # k=number of constraints=in_matrix_row_size+in_matrix_col_size

    # For both direct and iterative method
    p=array(0,dim=n)                                                                     # p vector allocation - input matrix linearization
    y=array(0,dim=n)                                                                     # y vector allocation - vector of "real" values
    c=array(0,dim=n)                                                                     # c vector allocation - reliability matrix linearization
    
    if (compress==0)
      V=array(0,dim=c(n,n))                                                              # V matrix allocation (Variance-CoVariance Matrix)
    else 
      V=Matrix(nrow = n, ncol = n, 0, sparse = TRUE)                                     # V matrix allocation (Variance-CoVariance Matrix) in optimized form
      
    if (compress==0)
      G=array(0,dim=c(k,n))                                                              # G matrix allocation (Constraints Matrix), k rows - n cols 
    else 
      G=Matrix(nrow = k, ncol = n, 0, sparse = TRUE)                                     # G matrix allocation (Variance-CoVariance Matrix) in optimized form

    h=array(0,dim=k)                                                                     # h vector
  
    # For the iterative method
    if (direct>0) {
      q=array(0,dim=k)                                                                     # q vector
      lambda_i=array(0.1,dim=k)                                                            # lambda vector   
      rho_i=array(0,dim=k)                                                                 # rho vector   
      phi_i=array(0,dim=k)                                                                 # phi vector   
      alpha_i=0                                                                            # initial alpha value
      beta_i=0                                                                             # initial beta value
      
      if (direct==2)
        if (compress==0)
          W=array(0,dim=c(k,k))                                                              # W matrix, pre-conditioning diagonal matrix     
        else
          W=Matrix(nrow = k, ncol = k, 0, sparse = TRUE)                                     # W matrix, pre-conditioning diagonal matrix in optimized form
    }
    
    out_matrix=array(0,dim=c(in_matrix_row_size, in_matrix_col_size))                      # allocates the output matrix for the balancing procedure
    out_row_marginals=t(array(0,dim=in_matrix_col_size))                                   # allocates the output column vector of row marginals
    out_col_marginals=array(0,dim=in_matrix_row_size)                                      # allocates the output row vector of col marginals
          
    diff_matrix=array(0,dim=c(in_matrix_row_size, in_matrix_col_size))                     # allocates the differences matrix in %
    diff_row_marginals=t(array(0,dim=in_matrix_col_size))                                  # allocates the differences column vector of row marginals in %
    diff_col_marginals=array(0,dim=in_matrix_row_size)                                     # allocates the differences row vector of col marginals in %

    discrepancy_row_marginals=t(array(0,dim=in_matrix_col_size))                           # difference between row sums and row marginals
    discrepancy_col_marginals=array(0,dim=in_matrix_row_size)                              # difference between row sums and col marginals
  
    # Definitions #
    if (simulation_test) {                                                                 # generates random matrixes
      in_matrix=matrix(runif(in_matrix_row_size*in_matrix_col_size)*max_random_matrix_element_values, ncol=in_matrix_col_size)

      # Generates a sparse matrix
      if ((type_of_simulation==2) || ((type_of_simulation=1) && (matrix_density<100)))    # generates sparse matrix
        for (i in 1:in_matrix_row_size) 
          for (j in 1:in_matrix_col_size) 
            if (runif(1)>(matrix_density/100))                                            # set to 0 all the values of the rest of the sparse matrix
              in_matrix[i,j]=0
      
      if (runif(1)<0.5) {
        row_marginals=rowSums(in_matrix)+matrix(runif(in_matrix_row_size)*max_random_marginals_perturbation, ncol=in_matrix_row_size)   # additional perturbation
      } else {
        row_marginals=rowSums(in_matrix)-matrix(runif(in_matrix_row_size)*max_random_marginals_perturbation, ncol=in_matrix_row_size)   # differential perturbation
      }

      if (runif(1)<0.5) {
        col_marginals=colSums(in_matrix)+matrix(runif(in_matrix_col_size)*max_random_marginals_perturbation, ncol=in_matrix_col_size)   # additional perturbation
      } else {
        col_marginals=colSums(in_matrix)-matrix(runif(in_matrix_col_size)*max_random_marginals_perturbation, ncol=in_matrix_col_size)   # differential perturbation
      }
    
    } else {
      if (load_from_xls) {                                                                # loaded from xls files
        in_matrix=xls_in_matrix
        row_marginals=xls_in_row_marginals
        col_marginals=xls_in_col_marginals
      } else {    
        # overloads the example matrices
        in_matrix=matrix(c(131.25,168.75,75,200,405,180,78.75,101.25,45), nrow = in_matrix_row_size, ncol = in_matrix_col_size, byrow = TRUE)
        row_marginals=c(400,850,250)
        col_marginals=c(500,700,350)
      }
    }

    if (simulation_test) {                                                               # generates random matrixes
      c_matrix=matrix(runif(in_matrix_row_size*in_matrix_col_size), ncol=in_matrix_col_size)
      c_row_marginals=matrix(runif(in_matrix_col_size), ncol=in_matrix_col_size)
      c_col_marginals=matrix(runif(in_matrix_row_size), ncol=in_matrix_row_size)
    } else {
      if (load_from_xls) {                                                                # loaded from xls files
          if (v_matrix) {   
            xls_tmp_in_matrix=xls_in_matrix                                               # create a temp matrix
            xls_tmp_in_row_marginals=xls_in_row_marginals                                 # create temp row marginals
            xls_tmp_in_col_marginals=xls_in_col_marginals                                 # create temp col marginals
            xls_tmp_in_matrix[xls_tmp_in_matrix==0]=1                                     # substitutes all the 0 values with 1 in order to avoid dividions by 0
            xls_tmp_in_row_marginals[xls_tmp_in_row_marginals==0]=1                       # substitutes all the 0 values with 1 in order to avoid dividions by 0
            xls_tmp_in_col_marginals[xls_tmp_in_col_marginals==0]=1                       # substitutes all the 0 values with 1 in order to avoid dividions by 0
            c_matrix=xls_v_matrix/xls_tmp_in_matrix                                       # make a right division
            c_row_marginals=xls_v_row_marginals/xls_tmp_in_row_marginals                  # make a right division
            c_col_marginals=xls_v_col_marginals/xls_tmp_in_col_marginals                  # make a right division
          } else {                                      
            c_matrix=xls_c_matrix
            c_row_marginals=xls_c_row_marginals
            c_col_marginals=xls_c_col_marginals
          }
      } else {    
        # overloads the example
        c_matrix=t(cbind(c(0.5,0.5,0.5),c(0.01,0.5,0.5),c(0.5,0.5,0.5)))                  # example of reliability matrix
        c_row_marginals=c(0.1,0.5,1.0)
        c_col_marginals=c(0.5,0.1,0.5)
      }
    }

    # Automatic building of G Matrix
    for (i in 1:in_matrix_row_size) 
      G[i,(1+((i-1)*in_matrix_col_size)):(in_matrix_col_size+((i-1)*in_matrix_col_size))]=rep(1,in_matrix_col_size)
    for (i in 1:in_matrix_row_size) 
      G[i,((in_matrix_col_size*in_matrix_row_size)+i)]=-1
    for (j in 1:in_matrix_col_size)
      for (i in 1:in_matrix_row_size)
        G[j+in_matrix_row_size,j+(i-1)*in_matrix_col_size]=1
    for (j in 1:in_matrix_col_size) 
      G[j+in_matrix_row_size,(in_matrix_row_size*in_matrix_col_size)+in_matrix_row_size+j]=-1
    
    ###################
    #   ELABORATION   #
    ###################

    # Pre-elaborations #
    # Linearization of input matrix, row marginals and column marginals and generation of p vector
    for (i in 1:in_matrix_row_size) 
      p[(1+((i-1)*in_matrix_col_size)):(in_matrix_col_size+((i-1)*in_matrix_col_size))]=in_matrix[i,]
    #this adds marginals    
    p[(1+(in_matrix_row_size*in_matrix_col_size)):(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size))]=row_marginals
    p[(1+(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size))):(in_matrix_col_size+(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size)))]=col_marginals

    # Linearization of reliability matrix (c_matrix), c_row_marginals and c_col_marginals and generation of c vector
    for (i in 1:in_matrix_row_size) 
      c[(1+((i-1)*in_matrix_col_size)):(in_matrix_col_size+((i-1)*in_matrix_col_size))]=c_matrix[i,]
    #this adds marginals    
    c[(1+(in_matrix_row_size*in_matrix_col_size)):(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size))]=c_row_marginals
    c[(1+(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size))):(in_matrix_col_size+(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size)))]=c_col_marginals

    #Load V matrix
    diag(V) = matrix(c * p)
                                                                                                   # %*% = matrix multiplication 
    # for the direct method
    p=matrix(p)                                                                                # p vector
    h=matrix(h)                                                                                # h vector
  
    # for the iterative method
    q=matrix(q)                                                                                # q vector
    lambda_i=matrix(lambda_i)                                                                  # lambda vector   
    rho_i=matrix(rho_i)                                                                        # rho vector   
    phi_i=matrix(phi_i)                                                                        # phi vector   

    #Stone's Estimator  (direct and iterative method)
    iter=0
 
    if (direct==0)                                                                             # Direct Method
      y=p-V %*% t(G) %*% solve(G %*% V %*% t(G)) %*% ((G %*% p) - h)                           # y=p-VG'(GVG')^-1(Gp-h)
    else {

      if (direct==1) {                                                                         # CG Method
        A=G %*% V %*% t(G)                                                                     # A symetric matrix, A = GVG'
        q=(G %*% p) - h                                                                        # q Vector, q=Gp-h
      } else if (direct==2) {                                                                  # CG with W matrix 
        # Pre-conditioning Matrix System
        diag(W)=matrix(1/sqrt(diag(G %*% V %*% t(G))))                                         # W matrix
        A=W %*% (G %*% V %*% t(G)) %*% W                                                       # A*=W(GVG')W 
        q=W %*% ((G %*% p) - h)                                                                # q*=Wq
      }
        
      # initial conditions i=0
      rho_i=q-A %*% lambda_i                                                                   # rho vector 
      phi_i=-rho_i                                                                             # phi i vector  (errore 1: phi era rho_i)
      alpha_i=((t(rho_i) %*% phi_i) / (t(phi_i) %*% A %*% phi_i))[1]                           # alpha_i value  (errore 2: c'era rho_i al posto del primo phi_i) 
      lambda_i=lambda_i + alpha_i * phi_i  
      beta_i=0                                                                                 # initial beta_i value (errore 5: beta non veniva inizializzato a 0      
    
      # main iteration loop
      for(iter in 1:k) {                                                                       # Conjugate Gradient Method Iteration, convergence within k cycles  
        rho_i=rho_i-alpha_i * A %*% phi_i                                                      # rho i+1
        print(paste(iter,": ",abs(colMeans(rho_i)),sep=""))
        if (abs(colMeans(rho_i))<exit_point)                                                   # stop condition 
          break;
        #if (abs(rho_i))<exit_point)                                                           # stop condition 
        #  break;
        beta_i=((t(rho_i) %*%  A %*% phi_i) / (t(phi_i) %*% A %*% phi_i))[1]                   # beta_i value  (errore 3: varie differenze in questa formula)      
        phi_i=-rho_i+beta_i * phi_i                                                            # phi i+1    (errore 4: qui il segno davanti a rho_i è negativo
        alpha_i=((t(rho_i) %*% phi_i) / (t(phi_i) %*% A %*% phi_i))[1]                         # alpha_i value  (errore 2: c'era rho_i al posto del primo phi_i) 
        lambda_i=lambda_i + alpha_i * phi_i  
      }
   
      if (direct==2) 
        lambda_i=W %*% lambda_i                                                                # re-adjusting lambda
        
      y=p-V %*% t(G) %*% lambda_i                                                              # y=p-VG'lambda, approximated Stone's formula
    }
                                                                                               # %*% = matrix multiplication 
    # Post elaborations #
    # De-Linearization of output matrix, row marginals and column marginals
    for (i in 1:in_matrix_row_size) 
      out_matrix[i,]=y[(1+((i-1)*in_matrix_col_size)):(in_matrix_col_size+((i-1)*in_matrix_col_size))]
    #this adds marginals    
    out_row_marginals=y[(1+(in_matrix_row_size*in_matrix_col_size)):(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size))]
    out_col_marginals=y[(1+(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size))):(in_matrix_col_size+(in_matrix_row_size+(in_matrix_row_size*in_matrix_col_size)))]

    # Calculating the differences matrix, row marginals and column marginals in %
    diff_matrix=(out_matrix-in_matrix)*100/out_matrix
    diff_row_marginals=(out_row_marginals-row_marginals)*100/out_row_marginals
    diff_col_marginals=(out_col_marginals-col_marginals)*100/out_col_marginals

    # Calculating the discrepancies between row and col sums of output matrix and marginals 
    discrepancy_row_marginals=rowSums(out_matrix)-out_row_marginals                            # difference between row sums and row marginals
    discrepancy_col_marginals=colSums(out_matrix)-out_col_marginals                            # difference between col sums and col marginals
  
    execution_time=proc.time() - ptm                                                           # Execution time
  
    if (type_of_simulation==1) {
      times_array[1,g,m]=execution_time[1]                                                     # capture the time in an array
      times_array[2,g,m]=as.integer(in_matrix_row_size)                                        # save the number of rows of matrix
      times_array[3,g,m]=as.integer(in_matrix_col_size)                                        # save the number of columns of matrix
      times_array[4,g,m]=as.integer(direct)                                                    # save the adopted method
      times_array[5,g,m]=iter                                                                  # save the number of iterations
    } else if (type_of_simulation==2) {
      times_array[1,g,m]=execution_time[1]                                                     # capture the time in an array
      times_array[2,g,m]=as.integer(matrix_density)                                            # save the matrix density
      times_array[4,g,m]=as.integer(direct)                                                    # save the adopted method
    }

    # Compute the error
    if (type_of_simulation==1) {
      errors_array[1,g,m]=colMeans(cbind(abs(discrepancy_row_marginals),0))[1]                 # Error of the process
      errors_array[2,g,m]=as.integer(in_matrix_row_size)                                       # save the number of rows of matrix
      errors_array[3,g,m]=as.integer(in_matrix_col_size)                                       # save the number of columns of matrix
      errors_array[4,g,m]=as.integer(direct)                                                   # save the adopted method
    } else if (type_of_simulation==2) {
      errors_array[1,g,m]=colMeans(cbind(abs(discrepancy_row_marginals),0))[1]                 # Error of the process
      errors_array[2,g,m]=as.integer(matrix_density)                                           # save the number of rows of matrix
      errors_array[4,g,m]=as.integer(direct)                                                   # save the adopted method
    }
  }     # end of g for
}    # end of m for

if (save_to_xls) {
  library(xlsx)
  print("SAVING OUTPUT MATRICES TO FILE .XLSX")
  res1=write.xlsx(out_matrix, paste(data_output,"/",experiment_name,"/balanced_matrix.xlsx", sep=""), sheetName="Balanced Output Matrix", col.names=FALSE, row.names=FALSE)
  res2=write.xlsx(out_row_marginals, paste(data_output,"/",experiment_name,"/balanced_out_row_marginals.xlsx", sep=""), sheetName="Balanced Output Row Marginals.xlsx", col.names=FALSE, row.names=FALSE)
  res3=write.xlsx(t(out_col_marginals), paste(data_output,"/",experiment_name,"/balanced_out_col_marginals.xlsx", sep=""), sheetName="Balanced Output Column Marginals.xlsx", col.names=FALSE, row.names=FALSE)
  res4=write.xlsx(in_matrix, paste(data_output,"/",experiment_name,"/unbalanced_matrix.xlsx", sep=""), sheetName="Unbalanced Output Matrix", col.names=FALSE, row.names=FALSE)
  res5=write.xlsx(t(row_marginals), paste(data_output,"/",experiment_name,"/unbalanced_out_row_marginals.xlsx", sep=""), sheetName="Unbalanced Output Row Marginals.xlsx", col.names=FALSE, row.names=FALSE)
  res6=write.xlsx(col_marginals, paste(data_output,"/",experiment_name,"/unbalanced_out_col_marginals.xlsx", sep=""), sheetName="Unbalanced Output Column Marginals.xlsx", col.names=FALSE, row.names=FALSE)
  res7=write.xlsx(c_matrix, paste(data_output,"/",experiment_name,"/reliability_matrix.xlsx", sep=""), sheetName="Reliability C Matrix", col.names=FALSE, row.names=FALSE)
  res8=write.xlsx(t(c_row_marginals), paste(data_output,"/",experiment_name,"/reliability_row_marginals.xlsx", sep=""), sheetName="Reliability C Row Marginals.xlsx", col.names=FALSE, row.names=FALSE)
  res9=write.xlsx(c_col_marginals, paste(data_output,"/",experiment_name,"/reliability_col_marginals.xlsx", sep=""), sheetName="Reliability C Column Marginals.xlsx", col.names=FALSE, row.names=FALSE)
}

#Plot of results and eventually save to files 
if (plot_all) {
	library(calibrate)										# Per usare textxy

	if (simulation_test>0)
    if (type_of_simulation==1) {
     
      ##### Times Plot #####
      
      # Determina il massimo dei dati sulle y
      max_y <- max(times_array[1,,])
    	max_y=max_y+max_y*lim_perc/100         							                                    # add a percentage to max_y	
      max_x=number_of_tests
   
      # Save plot on a file
      if (save_plot_to_file)                                                                                               
          dir.create(paste(data_output,"/",experiment_name,sep=""))                          	# Create the folder where saving all data			

      if (save_plot_to_file==1) 
        jpeg(paste(data_output,"/",experiment_name,"/stress_test_by_matrix_size.jpg", sep=""), width = 1280, height = 1024, bg="white", quality=100)	# Save all into jpeg     		
      if (save_plot_to_file==2)
  	 	  pdf(paste(data_output,"/",experiment_name,"/stress_test_by_matrix_size.pdf", sep=""), width = 15, height = 15)	                              # Save all into pdf
			
      if (sep_windows)
        dev.new()										               # Add a new window

      if (plot_type) {
        plot(times_array[1,,1], xlab="Input Matrix Size", ylab="Stone's Estimator Computation Times (secs)", xlim=c(0,max_x), ylim=c(0,max_y), type="b", lty=1, lwd=2, col=line_colors[1], axes=FALSE)		# Draw times line
    		#textxy(1:number_of_tests, times_array[1,,1],labs=times_array[5,,],m=c(0,0),cx=1.5)
        axis(1, las=1, at=0:number_of_tests, labels=c("0 x 0",paste(times_array[2,,1]," x ",times_array[3,,1])))
        axis(2, las=2, at=seq(0,max_y,by=max_y/number_of_tests), labels=round(seq(0,max_y,by=max_y/number_of_tests),2))
  
        # eventually adds other methods' data
        if (number_of_methods>1)
          for (m in 2:number_of_methods) 
            lines(times_array[1,,m], type="b", lty=1, lwd=2, col=line_colors[m], pch=sym_table[m])

        title(main = "Chart of Stone's Estimator Elaboration Times depending on Input Matrix Size Increases. Time is computated in seconds.\n", font.main = 4)
        legend(1, max_y, paste("Seconds of Time over Matrix Size, ",methods_names, sep=""), cex=1.6, col=line_colors, pch=sym_table, lty=1, lwd=2);
      } else {
        
        for (m in 1:number_of_methods) {
          times_tmp_array[m,]=times_array[1,,m]         # capture times of methods
          iters_tmp_array[m,]=times_array[5,,m]         # capture iterations of methods
        }       
        
      	barplot(times_tmp_array, xlab="Input Matrix Size", ylab="Stone's Estimator Computation Times (secs)", xlim=c(0,max_x), ylim=c(0,max_y), space=c(0.1,rep(rep(c(.1,2),c(number_of_methods-1,1)),number_of_tests-1),rep(0.1,number_of_methods-1)), beside=TRUE, names.arg=c(1:number_of_tests), col=line_colors)
        title(main = "Chart of Stone's Estimator Elaboration Times depending on Input Matrix Size Increases. Time is computated in seconds.\n", font.main = 4)
        legend(1, max_y, paste("Seconds of Time over Matrix Size, ",methods_names, sep=""), cex=1.6, col=line_colors, lty=1, lwd=2);
      }
       
      ##### Errors Plot #####

      # Determina il massimo dei dati sulle y
      max_y <- max(errors_array[1,,])
    	max_y=max_y+max_y*lim_perc/100         							                                    # add a percentage to max_y	
      max_x=number_of_tests
   
      # Save plot on a file
      if (save_plot_to_file==1) 
        jpeg(paste(data_output,"/",experiment_name,"/errors_by_matrix_size.jpg", sep=""), width = 1280, height = 1024, bg="white", quality=100)	# Save all into jpeg     		
      if (save_plot_to_file==2)
  	 	  pdf(paste(data_output,"/",experiment_name,"/errors_by_matrix_size.pdf", sep=""), width = 15, height = 15)	                              # Save all into pdf
			
      if (sep_windows)
        dev.new()										               # Add a new window

      plot(errors_array[1,,1], xlab="Input Matrix Size", ylab="Stone's Estimator Computation Error", xlim=c(0,max_x), ylim=c(0,max_y), type="b", lty=1, lwd=2, col=line_colors[1], axes=FALSE)		# Draw times line
      axis(1, las=1, at=0:number_of_tests, labels=c("0 x 0",paste(times_array[2,,1]," x ",times_array[3,,1])))
      axis(2, pos=0.6, las=2, at=seq(0,max_y,by=max_y/number_of_tests), labels=seq(0,max_y,by=max_y/number_of_tests))
  
      # eventually adds other methods' data
      if (number_of_methods>1)
        for (m in 2:number_of_methods) 
          lines(errors_array[1,,m], type="b", lty=1, lwd=2, col=line_colors[m], pch=sym_table[m])

      title(main = "Chart of Stone's Estimator Elaboration Error depending on Input Matrix Size Increases. \n", font.main = 4)
      legend(1, max_y, paste("Error over Matrix Size, ",methods_names, sep=""), cex=1.6, col=line_colors, pch=sym_table, lty=1, lwd=2);
      
    } else if (type_of_simulation==2) {
      ##### Times Plot #####
      
      # Determina il massimo dei dati sulle y
      max_y <- max(times_array[1,,])
    	max_y=max_y+max_y*lim_perc/100         							                                          # add a percentage to max_y	
      max_x=number_of_tests
   
      # Save plot on a file
      if (save_plot_to_file)                                                                                               
          dir.create(paste(data_output,"/",experiment_name,sep=""))                          	# Create the folder where saving all data			

      if (save_plot_to_file==1) 
        jpeg(paste(data_output,"/",experiment_name,"/stress_test_by_matrix_density.jpg", sep=""), width = 1280, height = 1024, bg="white", quality=100)	# Save all into jpeg     		
      if (save_plot_to_file==2)
  	 	  pdf(paste(data_output,"/",experiment_name,"/stress_test_by_matrix_density.pdf", sep=""), width = 15, height = 15)	                              # Save all into pdf
			
       if (sep_windows)
         dev.new()										               # Add a new window

      plot(times_array[1,,1], xlab="Input Matrix Density", ylab="Stone's Estimator Computation Times (secs)", xlim=c(0,max_x), ylim=c(0,max_y), type="b", lty=1, lwd=2, col=line_colors[1], axes=FALSE)		# Draw times line
      axis(1, las=1, at=0:number_of_tests, labels=c("0",times_array[2,,1]))
      axis(2, las=2, at=seq(0,max_y,by=max_y/number_of_tests), labels=round(seq(0,max_y,by=max_y/number_of_tests),2))

      # eventually adds other methods' data
      if (number_of_methods>1)
        for (m in 2:number_of_methods) 
          lines(times_array[1,,m], type="b", lty=1, lwd=2, col=line_colors[m], pch=sym_table[m])
  
      title(main = "Chart of Stone's Estimator Elaboration Times depending on Input Matrix Density Increases. Time is computated in seconds.\n", font.main = 4)
      legend(1, max_y, paste("Seconds of Time over Matrix Density, ",methods_names, sep=""), cex=1.6, col=line_colors, pch=sym_table, lty=1, lwd=2);

      
      ##### Errors Plot #####

      # Determina il massimo dei dati sulle y
      max_y <- max(errors_array[1,,])
    	max_y=max_y+max_y*lim_perc/100         							                                          # add a percentage to max_y	
      max_x=number_of_tests
   
      # Save plot on a file
      if (save_plot_to_file==1) 
        jpeg(paste(data_output,"/",experiment_name,"/errors_by_matrix_density.jpg", sep=""), width = 1280, height = 1024, bg="white", quality=100)	# Save all into jpeg     		
      if (save_plot_to_file==2)
  	 	  pdf(paste(data_output,"/",experiment_name,"/errors_by_matrix_density.pdf", sep=""), width = 15, height = 15)	                              # Save all into pdf
			
       if (sep_windows)
         dev.new()										               # Add a new window

      plot(errors_array[1,,1], xlab="Input Matrix Density", ylab="Stone's Estimator Computation Error", xlim=c(0,max_x), ylim=c(0,max_y), type="b", lty=1, lwd=2, col=line_colors[1], axes=FALSE)		# Draw times line
      axis(1, las=1, at=0:number_of_tests, labels=c("0",times_array[2,,1]))
      axis(2, pos=0.6, las=2, at=seq(0,max_y,by=max_y/number_of_tests), labels=seq(0,max_y,by=max_y/number_of_tests))

      # eventually adds other methods' data
      if (number_of_methods>1)
        for (m in 2:number_of_methods) 
          lines(errors_array[1,,m], type="b", lty=1, lwd=2, col=line_colors[m], pch=sym_table[m])
  
      title(main = "Chart of Stone's Estimator Elaboration Error depending on Input Matrix Density Increases. \n", font.main = 4)
      legend(1, max_y, paste("Error over Matrix Density, ",methods_names, sep=""), cex=1.6, col=line_colors, pch=sym_table, lty=1, lwd=2);
    }
}