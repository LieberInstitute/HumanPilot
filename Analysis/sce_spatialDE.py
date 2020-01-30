import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

import NaiveDE
import SpatialDE

import glob
samp_counts = glob.glob('/fastscratch/myscratch/shicks1/HumanPilot/sample_data/by_sample_id/*_counts.csv')
samp_meta = np.ravel([[x[:-10] + 'meta.csv'] for x in samp_counts])
samp_output = np.ravel([[x[:-10] + 'spatialDE_results.csv'] for x in samp_counts])
df = pd.DataFrame({'counts':samp_counts, 'meta':samp_meta, 'out':samp_output})

for index, row in df.iterrows():
    counts = pd.read_csv(row['counts'], index_col=0) # load counts
    counts = counts.T[counts.sum(axis=0) >= 5].T  # Filter practically unobserved genes
    sample_info = pd.read_csv(row['meta'], index_col=0) # load meta data with spatial coordinates
    norm_expr = NaiveDE.stabilize(counts.T).T # remove tech variation
    resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(sum)').T
    X = sample_info[['imagerow', 'imagecol']]
    now = datetime.now().strftime("%H:%M:%S")
    print(now)
    results = SpatialDE.run(X, resid_expr)
    now = datetime.now().strftime("%H:%M:%S")
    print(now)
    results.to_csv(row['out']) # Save spatial results
    print("Finished =", row['counts'])

