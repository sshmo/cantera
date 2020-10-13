import pandas as pd
import csv

dict_df = {}
with open("./w1/results.txt") as results:
    results_csv = csv.reader(results ,delimiter = ',')
    for row in results_csv:
        dict_df[row[0]] = row[1:]
df = pd.DataFrame(dict_df)

for column in df.columns:
    df[column] = df[column].astype('float64')


products = set(df.columns) - set(['Tsolid', 'h_conv', 'Tgas', 'Z'])

Mw_atoms = {"C":12, "O":16, "N":14, "H":1, "S":0, "A":0,"R":0}
MW ={}
for product in products:
    MW[product] = 0
    for i in range(len(product)):
        if product[i].isdigit() or product[i] in '()':
            continue
        if i < len(product)-1:
            if product[i+1].isdigit():
                MW[product] += int(product[i+1])*Mw_atoms[product[i]]
            else:
                MW[product] += Mw_atoms[product[i]]
        else:
            MW[product] += Mw_atoms[product[i]]
MW["AR"] = 40

products_molar = []
for product in products:
    products_molar.append(product + '_molar')
    df[product + '_molar'] = df[product]/MW[product]

df[products_molar] = df[products_molar].div(df[products_molar].sum(axis=1), axis=0)

main_products = ["H2_molar", "CO_molar","CH4_molar", "O2_molar","CO2_molar","H2O_molar"]
            
main_products_df = df[main_products]
main_products_df.to_csv('results.csv')





# import numpy as np

# df = pd.read_csv("./w1/results.txt", delimiter=",", index_col=None)

# use_cols = tuple([i for i in range(1, df.shape[1])])

# np_df = np.loadtxt("./w1/results.txt", delimiter=",", usecols=use_cols)

# columns = [df.columns[0]] + [i for i in df.iloc[:,0]]

# dft = pd.DataFrame(np_df).T

# dft.columns = columns

# dft.columns

# #------------------------------------