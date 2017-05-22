import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt

df1 = pd.read_csv('./_output/gauge00014.txt', skiprows=3, names=['dum1','t','h','hu','hv','eta'], sep=" ",skipinitialspace=True )
df2 = pd.read_csv('./_output_twin/gauge00014.txt', skiprows=3, names=['dum1','t','h','hu','hv','eta'], sep=" ",skipinitialspace=True )
print df1
fig,ax = plt.subplots(1,1)
df1.plot('t','eta',ax=ax)
df2.plot('t','eta',ax=ax)
plt.savefig('yup.png')
