#' ## PYTHON
# install.packages(reticulate)
#' library(reticulate)
#' Launch interactive python session from R Studio
# repl_python()

'''start python script'''

'''
My first Python script.
 
(this is a multi-line comment)
'''

# Enter 'exit' or 'quit' in the console to exit the REPL and return to R.

# library(...):
  #  import ...
 
# rnorm():
  #  random.normalvariate()

# Import the library
import random
# Library function
random.normalvariate()
# Assign values to the fuctions mu and sigma
random.normalvariate(mu=3, sigma=0.3)
# Assign to a variable
test0 = random.normalvariate(mu=3, sigma=0.3)
# Print the variables
test0
# Create a vector
range(2,8)
# Manually enter a python list
test02 = [2,3,4,5,6,7,8]
# Create a list by using brackets []
[ii for ii in range(2,8)]

x = test02
# Create a character value
test03 = 'x'
test03
print('x)
'x'

# Belows equivalence in R: for (ii in c(2:8)) {'x}
['x' for ii in range (2,8)]
# Belows equivalence in R: FOO <- c(); for(ii in 2:8) {FOO <- c(FOO,'x')}
for ii in range(2,8): ii
# Belows equivalence in R: rnorm(7,3,0.3)
[random.normalvariate(mu=3, sigma=0.3) for ii in range(2,8)]

''' October 6th, 2021 '''
# Install pandas
# Run py_install('pandas') in R console
import pandas as pd

#  To load a python libarary import FOO as BAR
dat0 = pd.read_csv('C:/Users/chapa/Documents/UTHSCSA/Course Work/FAll 2021/TSCI5230/Git Repository/over_expressed_genes.csv')
dat0
dat0.info()
dat0.head()

# Create a summary of the data
dat0.describe()

# Arrange the data by a column
dat0.sort_values('log2FoldChange')
dat0.sort_values(['log2FoldChange', 'padj'], ascending = [False,True])

# Subset a column
dat0.log2FoldChange
dat0['log2FoldChange']
# By column and row
dat0.padj
dat0.padj.between(3.839904e-28, 1.049713e-30, inclusive=True)
# loc is for names and logical expressions
dat0.loc[dat0.padj.between(1.049713e-30, 3.839904e-28, inclusive=True),'SYMBOL']
dat0.loc[dat0.padj.between(1.049713e-30, 3.839904e-28, inclusive=True),['SYMBOL','log2FoldChange']]
# iloc is for numeric selection of rows/columns ie. dat0.iloc[10:30,2:5]
dat0.iloc[10:30,2:5]
# Python indexs things starting from 0, so it skipped row 30 becuase 0 counted as a row, same with columns
dat0.iloc[10:30,1:5]

dat0['log2Fold Change'] = 42
dat0.head()
dat0['log2Fold Change']

# Operates on columns inside of a data frame
dat0.eval('pvalue+padj')

# For plotting py_install('matplotlib') in R terminal
import matplotlib.pyplot as plt
AP1=dat0.plot.scatter('log2FoldChange','padj')
plt.show()
AP1=dat0[dat0.SYMBOL=='FGR'].plot.scatter('log2FoldChange','padj', c='pink', label='FGR')
plt.show()
AP2=dat0[dat0.SYMBOL=='DPEP1'].plot.scatter('log2FoldChange','padj', c='blue', label='DPEP1', ax=AP1)
plt.show()

''' October 27th, 2021 '''
# Importing NumPy - allows for use of arrays in Python
import numpy as np
# Importing for linear regression modeling
from sklearn.linear_model import LinearRegression
b = np.array([1,5,7,8])
print(b)

# Create a matrix
A = np.array([[1,7,5,4], [7,9,4,5]])
print(A)

# Matrix with unmatched array sizes gets converted into a list
B = np.array([[1,7,5,4], [7,9,5]])

# Create a copy
C = A.copy()

# Other functions
A.shape
A.transpose()

# Arithmatics
b-3
A-b
A*b
A+5

# Fitting a linear regression model
# Use py$ in R to grab data from the python environment
# ex: py$iris <- r_to_py(iris)
# Side note - Vies(py) in R will show you the python environment
yy = np.array(iris["Petal.Length"])
yy
xx = np.array(iris[["Sepal.Length", "Sepal.Width"]])
# Assing linear regression function to an object
LRModel = LinearRegression()
LRModel.fit(xx, yy)
# When searching for things you can do with the object those in pink are not a function.
# For example....
LRModel.intercept_
LRModel.coef_
# To then remove the data from python in r use iris.py <- py$iris

# Combine an array by stacking
D = np.array([4,6,2,7])
E = np.vstack([C, D])
F = np.array([4,7,8])
# Switch F from rows to column
# Alternatively, F = np.array([4,7,8]).reshape(3,1)
# In R transpose is t()
F = F.reshape(3,1)
F
# Now you can stack the vectors vertically
G = np.hstack([E,F])
G
# Create a new vector that is 1 value more than every number in G
H = G+1
# Add -4 to every number in vector G
G += -4
G
