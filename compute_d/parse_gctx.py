from cmapPy.pandasGEXpress.parse import parse
import pandas as pd
import numpy as np
import os
os.chdir('TReS/data/LINCS2_beta/')

gctx_file = "level5_beta_trt_cp_n720216x12328.gctx"
print('start parsing')
gctx_data = parse(gctx_file)
print('parsing is ok')
# Access the expression matrix
exp_mat = gctx_data.data_df
exp_mat = pd.DataFrame(exp_mat)

# Access the sample information
pd.DataFrame(exp_mat.columns).to_csv('col.csv', index = 0, sep = '\t')
# Access the row (gene) information
pd.DataFrame(exp_mat.index).to_csv('row.csv', index = 0, sep = '\t')

exp_mat.to_csv('exp_mat.csv', index = 0, sep = ',')

print('ok')
