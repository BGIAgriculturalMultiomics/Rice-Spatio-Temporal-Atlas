from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles,
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd

A=pd.read_csv("At_to_ZH.txt",sep='\t',index_col=0,header=None)
B=pd.read_csv("ZH_to_At.txt",sep='\t',index_col=0,header=None)

fn1="Rice.S01-S05_intergration.h5ad"
fn2="At.rossete_combine_pr.h5ad"

sam1=SAM()
sam1.load_data(fn1)

sam2=SAM()
sam2.load_data(fn2)

sams = {'ZH':sam1,'At':sam2}
sm = SAMAP(
        sams,
        f_maps="./SC/rice_leaf/maps/",
        )

sm.run(pairwise=True)
samap = sm.samap

sm = SAMAP(
        sams,
        f_maps = './SC/rice_leaf/maps/',
        resolutions = {'ZH':1,'At':1},
        save_processed=True
    )
sm.run(pairwise=True)
samap = sm.samap
samap.save_anndata("RiceAt.mg.h5ad")

sm.sams['ZH'].adata.obs['sample'].head()
sm.sams['At'].adata.obs['stage'].head()

sm = SAMAP(
        sams,
        f_maps = './SC/rice_leaf/maps/',
        keys = {'ZH':'sample','At':'stage'},
        resolutions = {'ZH':1,'At':1},
        save_processed=True
        )
sm.run(pairwise=True)
samap = sm.samap

keys = {'ZH':'sample','At':'stage'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
D.to_csv("RiceAt.mapping.scores.csv",index=False)
MappingTable.to_csv('RiceAt.MappingTable.new.csv',index=False)