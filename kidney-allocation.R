library(gurobi)
library(igraph)

############################### DATA PREPARATION ######################################

set.seed(0017)

# Generate number of patient-donor pairs
N = sample(20:25, 1)
N

# Generate compatibility matrix A1
S = c(1:N)
A1 = matrix(, nrow = 0, ncol = N)
for(i in 1:N) {
  A1 = rbind(A1, sample(S))
}
A1

# Generate compatibility matrix A2
A2 = matrix(, nrow = 0, ncol = N)
for(i in 1:N) {
  A2 = rbind(A2, round(runif(N, 0, 1)))
}
A2

################################### PART 1 ##########################################

# ------------------------- MINIMIZING TOTAL PREFEERNCE -----------------------------

# Matrix A1 has preferences where Pij = preference of recipient i to receive kidney from donor j 
# Xij = 1 if recipient i receives kidney of donor j, 0 otherwise
# Decision variables are in the following order: x11, x12, x13, ... x1n, x21, x22, ... x2n, ... xn1, xn2, ... xnn
# Assumption: A patient can receive from their own donor

# Set the objective function vector
obj = as.vector(t(A1))

# Set the variable types. There are 2 options: 1) Use binary variables 2) Use constraint xij>=0
vtype = matrix('B', nrow = 1, ncol = N*N)

# Set the operators vector
operators = matrix('=', nrow = 2*N, ncol = 1 )

# Set the A matrix 
A = matrix(nrow = 0, ncol = N*N)

# First N constraints are for: each recipient receives exactly one kidney
for(i in 1:N) {
  rowVector = vector("numeric", N*N)
  for(j in 1:N) {
    rowVector[j + (i-1)*N] = 1
  }
  A = rbind(A, rowVector)
}

# Last N constraints are for: each kidney is donated to exactly one recipient
for(i in 1:N) {
  rowVector = vector("numeric", N*N)
  for(j in 1:N) {
    rowVector[i + (j-1)*N] = 1
  }
  A = rbind(A, rowVector)
}

# Set the B vector
b = matrix(1, nrow = 2*N, ncol = 1 )

# Solve
model = list()
model$A = A
model$obj = obj
model$modelsense = "min"
model$rhs = b
model$sense = operators
model$vtype = vtype
result = gurobi(model)

# Get the values of the decision variables
solutionVarsMatrix = matrix(result$x, nrow = N,ncol = N, byrow = TRUE)
solutionVarsMatrix

preferences = vector("numeric",0)
for(i in 1:N) {
  rowVar = solutionVarsMatrix[i,]
  rowA1 = A1[i,]
  pref = sum(rowVar*rowA1)
  preferences = c(preferences, pref)
}
# Preferences for all nodes
preferences
# Total preference
sum(preferences)
# Max preference used in solution
max(preferences)

# Plot the graph
network = graph_from_adjacency_matrix(t(solutionVarsMatrix), mode="directed")
png("part1-minTotalPref.png", width=1200, height=1200)
par(mar=c(2,2,2,2))
plot.igraph(network)
title("Kidney Exchange Solution for Minimizing Total Preference",cex.main=2)
dev.off()

# ------------------------ MINIMIZING THE MAXIMUM PREFEERNCE --------------------------------

# Matrix A1 has preferences where Pij = preference of recipient i to receive kidney from donor j 
# Xij = 1 if recipient i receives kidney or donor j, 0 otherwise
# Decision variables are in the following order: x11, x12, x13, ... x1n, x21, x22, ... x2n, ... xn1, xn2, ... xnn, z
# where z is the maximum preference used in the solution

# Set the objective function vector
obj = matrix(0, nrow = 1, ncol = N*N + 1)
obj[N*N + 1] = 1

# Set the variable types. There are 2 options: 1) Use binary variables 2) Use constraint xij>=0
vtype = matrix('B', nrow = 1, ncol = N*N)
vtype = cbind(vtype, matrix('C', nrow=1, ncol=1))

# Set the operators vector
operators = matrix('=', nrow = 2*N, ncol = 1 )
operators = rbind(operators, matrix(">=", nrow=N, ncol = 1))

# Set the A matrix 
A = matrix(nrow = 0, ncol = N*N + 1)

# First N constraints are for: each recipient receives exactly one kidney
for(i in 1:N) {
  rowVector = vector("numeric", N*N + 1)
  for(j in 1:N) {
    rowVector[j + (i-1)*N] = 1
  }
  A = rbind(A, rowVector)
}

# Middle N constraints are for: each kidney is donated to exactly one recipient
for(i in 1:N) {
  rowVector = vector("numeric", N*N + 1)
  for(j in 1:N) {
    rowVector[i + (j-1)*N] = 1
  }
  A = rbind(A, rowVector)
}

# Last N constraints for for getting the maximum of the preferences used in the solution
for(i in 1:N) {
  rowVector = vector("numeric", N*N + 1)
  rowVector[N*N + 1] = 1
  for(j in 1:N) {
    rowVector[j + (i-1)*N] = -1*A1[i,j]
  }
  A = rbind(A, rowVector)
}
# write.csv(A, file="A.csv")

# Set the B vector
b = matrix(1, nrow = 2*N, ncol = 1)
b = rbind(b, matrix(0, nrow = N, ncol = 1))

# Solve
model = list()
model$A = A
model$obj = obj
model$modelsense = "min"
model$rhs = b
model$sense = operators
model$vtype = vtype
result = gurobi(model)

solutionVarsMatrix = matrix(head(result$x, -1), nrow = N, ncol = N, byrow = TRUE)
solutionVarsMatrix

preferences = vector("numeric",0)
for(i in 1:N) {
  rowVar = solutionVarsMatrix[i,]
  rowA1 = A1[i,]
  pref = sum(rowVar*rowA1)
  preferences = c(preferences, pref)
}
preferences
# Total preferences
sum(preferences)
# Max preference used in solution
max(preferences)

# Plot the graph
network = graph_from_adjacency_matrix(t(solutionVarsMatrix), mode="directed" )
png("part1-minMaxPref.png", width=1200, height=1200)
par(mar=c(2,2,2,2))
plot.igraph(network)
title("Kidney Exchange Solution for Minimizing Maximum Preference",cex.main=2)
dev.off()


# ------------------------------- 2-WAY EXCHANGES --------------------------------------

# Matrix A2 represents compatibility between patient i and donor j
# Make a graph where 1 represents a compatible arc (where patient i is compatible with donor j AND patient j is compatible with donor i)
# This graph is symmetric, but we only really need the upper triangle
graph = matrix(0, nrow=N, ncol=N)
for(i in 1:N) {
  for(j in 1:N) {
    # if i can receive from j and j can receive from i
    if(A2[i,j] == 1 && A2[j,i] == 1) {
      graph[i,j] =1
    }
  }
}

# We only need the upper triangle of the graph since we do not have a directed graph anymore and the graph is symmetric

# Decision variables are in the following order: x12, x13, x14... x1n, x23, x24 ... x2n, x34... x3n,... x(n-1)n
numVars = 0
for(i in 1:N-1) {
  numVars = numVars + i
}

# Set the objective function vector
# Use a decision variable for each possible pair of nodes for simplicity of coding. We will set decision vars to zero 
# for pairs of nodes that are not two-way compatible
obj = matrix(1, nrow = 1, ncol = numVars)

# Set the variable types
vtype = matrix('B', nrow = 1, ncol = numVars)

# Set the operators vector
operators = matrix('<=', nrow = N, ncol = 1)

# Set the B vector
b = matrix(1, nrow = N, ncol = 1 )

# Set the A matrix 
A = matrix(nrow = 0, ncol = numVars)

# Add constraints for each node only being involved in one exchange (one constraint per row/node)
for(i in 1:N) {
  temp = matrix(0, nrow=N, ncol=N)
  for(j in i:N) {
    temp[i,j] = 1
  }
  for(k in 1:(i-1)) {
    temp[k,i] = 1
  }
  temp[i,i]=0
  # temp matrix will have an L-shaped pattern that represents the variables we want to include in our constraint
  # our decision variables follow the pattern of the upper-diagonal matrix so we can use our temp matrix to make our constraints
  rowVector = vector("numeric", 0)
  for(l in 1:(N-1)) {
    tempRow = as.vector(temp[l,])
    rowVector = c(rowVector, tempRow[(l+1):N])
  }
  A = rbind(A, rowVector)
}

# ADD CONSTRAINTS TO MAKE XIJ=0 IF NOT IN GRAPH
# Iterate through the upper triangle of the graph not including the diagonal
# Keep track of the index of the corresponding decision variable
indexOfDecisionVar = 0
# Iterate through rows
for(i in 1:(N-1)) {
  # Iterate through columns starting with the element one to the right in the diagonal
  for(j in (i+1):N) {
    indexOfDecisionVar = indexOfDecisionVar + 1
    if(graph[i,j] == 0) {
      # Add to A
      rowVector = vector("numeric", numVars)
      rowVector[indexOfDecisionVar] = 1
      A = rbind(A, rowVector)
      # Add to operators
      operators = rbind(operators, matrix('=', nrow = 1, ncol = 1))
      # Add to b
      b = rbind(b, matrix(0, nrow=1, ncol=1))
    }
  }
}

# Solve
model = list()
model$A = A
model$obj = obj
model$modelsense = "max"
model$rhs = b
model$sense = operators
model$vtype = vtype
result = gurobi(model)

resultVars = result$x
resultMatrix = matrix(0,nrow=0,ncol=N)
for(i in N:2) {
  zeroes = vector("numeric", N-i+1)
  resultMatrix = rbind(resultMatrix, c(zeroes, resultVars[1:(i-1)]))
  resultVars = resultVars[i:length(resultVars)]
}
lastRow = matrix(0, nrow = 1, ncol = N)
resultMatrix = rbind(resultMatrix, lastRow)
resultMatrix

# Plot the graph
network = graph_from_adjacency_matrix(resultMatrix, mode="upper" )
png("part1-twoWay.png", width=1200, height=1200)
par(mar=c(2,2,2,2))
plot.igraph(network)
title("Kidney Exchange Solution for Two-Way Exchanges", cex.main=2)
dev.off()


################################### PART 2 ##########################################

# Matrix A1 represents preferences between patient i and donor j. Preference of 10 or less means compatible.
# Make a graph where 1 represents a compatible arc (where patient i is compatible with donor j AND patient j is compatible with donor i)
# This graph is symmetric, but we only really need the upper triangle
graph = matrix(0, nrow=N, ncol=N)
for(i in 1:N) {
  for(j in 1:N) {
    # if i can receive from j and j can receive from i (both preferences are 10 or under)
    if(A1[i,j] <= 10 && A1[j,i] <=10) {
      graph[i,j] =1
    }
  }
}

# We only need the upper triangle of the graph since we do not have a directed graph anymore and the graph is symmetric

# Decision variables are in the following order: x12, x13, x14... x1n, x23, x24 ... x2n, x34... x3n,... x(n-1)n
numVars = 0
for(i in 1:N-1) {
  numVars = numVars + i
}

# Set the objective function vector
# Use a decision variable for each possible pair of nodes for simplicity of coding. We will set decision vars to zero 
# for pairs of nodes that are not two-way compatible
obj = matrix(1, nrow = 1, ncol = numVars)

# Set the variable types
vtype = matrix('B', nrow = 1, ncol = numVars)

# Set the operators vector
operators = matrix('<=', nrow = N, ncol = 1)

# Set the B vector
b = matrix(1, nrow = N, ncol = 1 )

# Set the A matrix 
A = matrix(nrow = 0, ncol = numVars)


# Add constraints for each node only being involved in one exchange (one constraint per row/node)
for(i in 1:N) {
  temp = matrix(0, nrow=N, ncol=N)
  for(j in i:N) {
    temp[i,j] = 1
  }
  for(k in 1:(i-1)) {
    temp[k,i] = 1
  }
  temp[i,i]=0
  # temp matrix will have an L-shaped pattern that represents the variables we want to include in our constraint
  # our decision variables follow the pattern of the upper-diagonal matrix so we can use our temp matrix to make our constraints
  rowVector = vector("numeric", 0)
  for(l in 1:(N-1)) {
    tempRow = as.vector(temp[l,])
    rowVector = c(rowVector, tempRow[(l+1):N])
  }
  A = rbind(A, rowVector)
}

# ADD CONSTRAINTS TO MAKE XIJ=0 IF NOT IN GRAPH
# Iterate through the upper triangle of the graph not including the diagonal
# Keep track of the index of the corresponding decision variable
indexOfDecisionVar = 0
# Iterate through rows
for(i in 1:(N-1)) {
  # Iterate through columns starting with the element one to the right in the diagonal
  for(j in (i+1):N) {
    indexOfDecisionVar = indexOfDecisionVar + 1
    if(graph[i,j] == 0) {
      # Add to A
      rowVector = vector("numeric", numVars)
      rowVector[indexOfDecisionVar] = 1
      A = rbind(A, rowVector)
      # Add to operators
      operators = rbind(operators, matrix('=', nrow = 1, ncol = 1))
      # Add to b
      b = rbind(b, matrix(0, nrow=1, ncol=1))
    }
  }
}

# Solve
model = list()
model$A = A
model$obj = obj
model$modelsense = "max"
model$rhs = b
model$sense = operators
model$vtype = vtype
result = gurobi(model)

resultVars = result$x
resultMatrix = matrix(0,nrow=0,ncol=N)
for(i in N:2) {
  zeroes = vector("numeric", N-i+1)
  resultMatrix = rbind(resultMatrix, c(zeroes, resultVars[1:(i-1)]))
  resultVars = resultVars[i:length(resultVars)]
}
resultMatrix

optimalCardinality = result$objval
optimalCardinality

# Start adding constraints to see if we can prioritize patients while remaining optimal

# Funcion to add constraint based on the index of the priority patient
addConstraint = function(priorityPatient){
  # Add constraint for priority patient
  temp = matrix(0, nrow=N, ncol=N)
  for(j in priorityPatient:N) {
    temp[priorityPatient,j] = 1
  }
  for(k in 1:(priorityPatient-1)) {
    temp[k,priorityPatient] = 1
  }
  temp[priorityPatient,priorityPatient] = 0
  # temp matrix will have an L-shaped pattern that represents the variables we want to include in our constraint
  # our decision variables follow the pattern of the upper-diagonal matrix so we can use our temp matrix to make our constraints
  rowVector = vector("numeric", 0)
  for(l in 1:(N-1)) {
    tempRow = as.vector(temp[l,])
    rowVector = c(rowVector, tempRow[(l+1):N])
  }
  # return the constraint to be added
  return(rowVector)
}
  
solveModel = function(){
  # Solve model
  model = list()
  model$A = A
  model$obj = obj
  model$modelsense = "max"
  model$rhs = b
  model$sense = operators
  model$vtype = vtype
  result = gurobi(model)
  
  # Print the results if the solution is optimal (don't print if solution is infeasible)
  if(result$status == "OPTIMAL") {
    resultVars = result$x
    resultMatrix = matrix(0,nrow=0,ncol=N)
    for(i in N:2) {
      zeroes = vector("numeric", N-i+1)
      resultMatrix = rbind(resultMatrix, c(zeroes, resultVars[1:(i-1)]))
      resultVars = resultVars[i:length(resultVars)]
    }
    print("Result Matrix:")
    print(resultMatrix)
    
    # Get all the preferences for this solution
    preferences = vector("numeric",0)
    for(i in 1:(N-1)) {
      for(j in 1:N) {
        if(resultMatrix[i,j] == 1) {
          prefA = A1[i,j]
          prefB = A1[j,i]
          preferences = c(preferences, prefA, prefB)
        }
      }
    }
    # Preferences are grouped by pairs, not ordered by patient number 
    print("Preferences:")
    print(preferences)
    print("Total Preference:")
    print(sum(preferences))
    print("Max Preference:")
    print(max(preferences))
  }
  return(result)
}

# Iteratively add constraints based on priority
for(i in 1:N) {
  # Add constraint and solve
  constraint = addConstraint(i)
  A = rbind(A, constraint)
  operators = rbind(operators, matrix("=", nrow=1, ncol=1))
  b = rbind(b, matrix(1, nrow=1, ncol=1))
  modelResult = solveModel()
  # cardinality will be NULL if model is infeasible
  cardinality = modelResult$objval
  status = modelResult$status
  # remove constraint if constraint causes infeasibility or a cardinality less than in the optimal solution
  if(status == "INFEASIBLE" || status == "INF_OR_UNBD" || cardinality < optimalCardinality) {
    # remove constraint from A
    A = A[1:nrow(A)-1,]
    # remove from b
    b = head(b, -1)
    # remove from operators 
    operators = head(operators, -1)
  } 
}

# Get the final solution
modelResult = solveModel()

# Build the graph
resultVars = modelResult$x
resultMatrix = matrix(0,nrow=0,ncol=N)
for(i in N:2) {
  zeroes = vector("numeric", N-i+1)
  resultMatrix = rbind(resultMatrix, c(zeroes, resultVars[1:(i-1)]))
  resultVars = resultVars[i:length(resultVars)]
}
lastRow = matrix(0, nrow = 1, ncol = N)
resultMatrix = rbind(resultMatrix, lastRow)
network = graph_from_adjacency_matrix(resultMatrix, mode="upper" )
png("part2-twoWay.png", width=1200, height=1200)
par(mar=c(2,2,2,2))
plot.igraph(network)
title("Kidney Exchange Solution for Two-Way Exchanges with the Trading Algorithm", cex.main=2)
dev.off()
