##
##  satAttack.r
##
##  Demonstrate SAT solver on table of aggregate releases for reconstruction attack
##
##  jH 2019.2.4
##


# Here is our simple dataset:
#  V1    V2
#  1: 1  6:0
#  2: 1  7:1
#  3: 0  8:1
#  4: 0  9:1
#  5: 0  10:0


V1 <- c(1,1,0,0,0)
V2 <- c(0,1,1,1,0)
data <- data.frame(V1, V2)
print(data)


# solve a boolean formula
    # each variable is an integer
    # negations are negative integers
    # each vector is a disjunction
    # the list is a conjunction
    formula <- list(
      ## Sum V1 = 2
      # Of any four values of V1, at least one must be TRUE
      c(1L, 2L, 3L, 4L), c(1L, 2L, 3L, 5L), c(1L, 2L, 4L, 5L), c(1L, 3L, 4L, 5L), c(2L, 3L, 4L, 5L), 
      # Of any four values of V1, at least one must be FALSE
      c(-1L, -2L, -3L), c(-1L, -2L, -4L), c(-1L, -2L, -5L), c(-1L, -3L, -4L), c(-1L, -3L, -5L), c(-1L, -4L, -5L), c(-2L, -3L, -4L), c(-2L, -3L, -5L), c(-2L, -4L, -5L), c(-3L, -4L, -5L), 
      ## Sum V2 = 3
      # Of any three values of V2, at least one must be TRUE
      c(6L, 7L, 8L), c(6L, 7L, 9L), c(6L, 7L, 10L), c(6L, 8L, 9L), c(6L, 8L, 10L), c(6L, 9L, 10L), c(7L, 8L, 9L), c(7L, 8L, 10L), c(7L, 9L, 10L), c(8L, 9L, 10L),
      # Of any four values of V2, at least one must be FALSE
      c(-6L, -7L, -8L, -9L), c(-6L, -7L, -8L, -10L), c(-6L, -7L, -9L, -10L), c(-6L, -8L, -9L, -10L), c(-7L, -8L, -9L, -10L) 
      )


    # find the constraints of n objects which sum to s
    sumConstraint <- function(low, high, sum){
      constraint1 <- combn(low:high, (high-low+1)-sum+1, simplify=FALSE)
      constraint2 <- combn(-low:-high, sum+1, simplify=FALSE)
      constraints <- c(constraint1, constraint2)
      return(constraints)
    }

    formula2 <- c(sumConstraint(low=1, high=5, sum=2)), sumConstraint(low=6, high=10, sum=3))

    cat("check formulas contain same terms \n")
    print(formula %in% formula2)
    print(formula2 %in% formula)


    # find the formula for a dataset that defines the variable sums
    generateFormula <- function(data){
      formula <- NULL
      n <- nrow(data)
      for(i in 1:ncol(data)){
        formula <- c(formula, sumConstraint(low=((i-1)*n+1), high=i*n, sum=sum(data[,i])))
      }
      return(formula)
    }

    formula3 <- generateFormula(data)

    cat("check formulas contain same terms \n")   
    print(identical(formula,formula3))


    res <- picosat_sat(formula)
    # If we need to set a variable to a fixed value
    # e.g. 1L = TRUE and 2L = FALSE
    #res <- picosat_sat(formula, assumptions = c(1L, -2L))
    
    picosat_solution_status(res)
    
    # get further information about the solution process
    picosat_variables(res)
    picosat_added_original_clauses(res)
    picosat_decisions(res)
    picosat_propagations(res)
    picosat_visits (res)
    picosat_seconds(res)

