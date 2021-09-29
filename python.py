#' ## PYTHON
# install.packages(reticulate)
# library(reticulate)
#' Launch interactive python session from R Studio
# repl_python()

'''start python script'''

'''
My first Python script.
 
(this is a multi-line comment)
'''

'Enter 'exit' or 'quit' in the console to exit the REPL and return to R.

'library(...):
  '  import ...

'rnorm():
  '  random.normalvariate()

' Import the library
import random
' Library function
random.normalvariate()
' Assign values to the fuctions mu and sigma
random.normalvariate(mu=3, sigma=0.3)
' Assign to a variable
test0 = random.normalvariate(mu=3, sigma=0.3)
' Print the variables
test0
'Create a vector
range(2,8)
'Manually enter a python list
test02 = [2,3,4,5,6,7,8]
' Create a list by using brackets []
[ii for ii in range(2,8)]

x = test02
'Create a character value
test03 = 'x'
test03
print('x)
'x'

' Belows equivalence in R: for (ii in c(2:8)) {'x}
['x' for ii in range (2,8)]
' Belows equivalence in R: FOO <- c(); for(ii in 2:8) {FOO <- c(FOO,'x')}
for ii in range(2,8): ii
' Belows equivalence in R: rnorm(7,3,0.3)
[random.normalvariate(mu=3, sigma=0.3) for ii in range(2,8)]
