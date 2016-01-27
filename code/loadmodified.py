import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
import os


def load_1000_funct_connectome(data_folder='data', location='Baltimore',
                               verbose=True):
    """Load the functional connectivity dataset 1000_functional
    connectomes available at:
    http://umcd.humanconnectomeproject.org/umcd/default/browse_studies
    Parameters:
    ----------
    data_folder: string
                Path to the folder containing all the data files
    verbose: bool
    Return:
    ------
    X, y: (ndarray, array)
         The dataset with the class labels
    """
    path_folder = os.path.join(data_folder, 'Functional_Connectomes',
                               'Locations', location)
    desc_file = os.path.join(data_folder, 'Functional_Connectomes',
                             '1000_Functional_Connectomes.csv')

    desc = pd.read_csv(desc_file, sep=',')
    name = desc['upload_data.network_name']
    pool = desc['upload_data.gender']

    dt_name = dict()
    for i, v in enumerate(name):
        dt_name[v] = i

    dt_cls = dict()
    cs = np.unique(pool)
    for i, v in enumerate(cs):
        dt_cls[v] = i

    dirs = os.listdir(path_folder)
    dirs.sort()
    connects = []

    names = []
    y = []

    for i, v in enumerate(dirs):
        spl = v.split('_')
        if spl[-2] == 'matrix':
            connects.append(np.loadtxt(os.path.join(path_folder, v)))
            nm = '_'.join(spl[:-3])
            names.append(nm)
            y.append(dt_cls[pool[dt_name[nm]]])
            assert cs[y[-1]] == pool[dt_name[nm]], 'wrong class assignment'

    X = np.zeros((len(connects), (connects[0].shape[0] *
                                  (connects[0].shape[1] - 1)) / 2))
    triu_idx = np.triu_indices(connects[0].shape[0], k=1)
    for i, v in enumerate(connects):
        X[i] = v[triu_idx]

    # Dealing with negative correlations
    X = np.where(X > 0, X, 0.)
    # X = np.abs(X)

    # Being sure that labels are 0-1
    le = LabelEncoder()
    y = le.fit_transform(y)

    # Sorting data by labels
    ast = np.argsort(y)
    X = X[ast]
    y = y[ast]

    mm = MinMaxScaler()
    X = mm.fit_transform(X.T).T

    if verbose:
        print 'n_samples: %s, n_samples_by_class: (%s - %s)' % (len(y),
                                                                len(y[y == 0]),
                                                                len(y[y == 1]))

    return X, y
