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

# Belows equivalence in R: for (ii in c(2:8)) {'x}
['x' for ii in range (2,8)]
# Belows equivalence in R: FOO <- c(); for(ii in 2:8) {FOO <- c(FOO,'x')}
for ii in range(2,8): ii
# Belows equivalence in R: rnorm(7,3,0.3)
[random.normalvariate(mu=3, sigma=0.3) for ii in range(2,8)]

''' October 6th, 2021 '''
' Install pandas
# Run py_install('pandas') in R console
import pandas as pd

'  To load a python libarary import FOO as BAR
dat0 = pd.read_csv('C:/Users/chapa/Documents/UTHSCSA/Course Work/FAll 2021/TSCI5230/Git Repository/over_expressed_genes.csv')
dat0
dat0.info()
dat0.head()

' Create a summary of the data
dat0.describe()

' Arrange the data by a column
dat0.sort_values('log2FoldChange')
dat0.sort_values(['log2FoldChange', 'padj'], ascending = [False,True])

' Subset a column
dat0.log2FoldChange
dat0['log2FoldChange']
' By column and row
dat0.padj
dat0.padj.between(3.839904e-28, 1.049713e-30, inclusive=True)
' loc is for names and logical expressions
dat0.loc[dat0.padj.between(1.049713e-30, 3.839904e-28, inclusive=True),'SYMBOL']
dat0.loc[dat0.padj.between(1.049713e-30, 3.839904e-28, inclusive=True),['SYMBOL','log2FoldChange']]
' iloc is for numeric selection of rows/columns ie. dat0.iloc[10:30,2:5]
dat0.iloc[10:30,2:5]
' Python indexs things starting from 0, so it skipped row 30 becuase 0 counted as a row, same with columns
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
