import os

oak_base = globals().get('oak_base', os.environ["OAK_BASE"])
base_path = globals().get('base_path', os.path.join(oak_base,'users/daphna/eranb'))
phenotypes_df = globals().get('phenotypes_df',os.path.join(base_path,'Eran_Bmal1_Meta data analysis_05.18.2020_corrected.xlsx'))
analyses_path = globals().get('analyses_path',os.path.join(base_path,'analyses'))
y_category = globals().get('y_category','Genotype')
boolean_types = globals().get('boolean_types',{'Genotype':{'KO':1,'WT':0},'Sex':{'Male':0,'Female':1}})
