inputCond <- function(init = 0, prompt = "Please enter a value for x: ", variable = "x", condition = "x < 3"){
  #
  # Takes user input and tests it against a given condition.
  #
  # init      Numeric. The initial value of the variable
  # prompt    Character string. Printed to prompt the user for input
  # variable  Character string. The name of the variable being tested in the condition
  # condition Character string. The logical test of the variable. If true, the user input is returned.
  #
  x <- variable
  assign(x, init)
  while (x == init){
    
    x <- as.numeric(readline(paste(prompt)))
    if (eval(parse(text = condition)) == F){
      cat(paste("Please choose another value such that ", condition, sep = ""))
      x <- init
    }
    else {
      return(x)
    }
  }
}
