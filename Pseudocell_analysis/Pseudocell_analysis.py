import numpy as np
import pandas as pd

def Pseudocell_analysis_pipeline(DGE_tab_data, phenotype_tab_data, pseudocell_size=20, discard_t=0.8):
    """
    @inputs
    DGE_tab_data: cells * genes matrix
    phenotype_tab_data: should contain ordered CellID, Celltype
    pseudocell_size: size of pseudocell                                 DEFAULT: 20
    discard_t: cell number < discard_t * pseudocell will be discard     DEFAULT: 0.8
    """

    data = pd.read_csv(DGE_tab_data, sep=',', index_col=0, header=0).T
    anno = pd.read_csv(phenotype_tab_data, sep =',', index_col=None, header=0)

    anno.columns = ["CellID", "Tissue", "Celltype"]
    anno.index = anno["CellID"].values

    anno['pseudo.id'] = 0

    for cell_type in np.unique(anno['Celltype']):
        idx = anno["Celltype"]==cell_type
        anno.loc[idx, "pseudo.id"] = range(np.sum(idx))

    anno['pseudo.id'] = np.floor_divide(anno['pseudo.id'], pseudocell_size).astype(str)
    anno['pseudo.id'] = anno["Celltype"] +'_Cell'+ anno['pseudo.id']
    idx = anno.groupby(['pseudo.id'])['CellID'].count() >= discard_t * pseudocell_size

    data['pseudo.id'] = anno['pseudo.id'].values
    data_mean = data.groupby(['pseudo.id']).mean()[idx]

    return data_mean