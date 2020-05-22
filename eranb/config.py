import os

oak_base = globals().get('oak_base', os.environ["OAK_BASE"])
base_path = globals().get('base_path', os.path.join(oak_base,'users/daphna/eranb'))
phenotypes_df = globals().get('phenotypes_df',os.path.join(base_path,'Eran_Bmal1_Meta data analysis_05.18.2020_corrected 2.xlsx'))
analyses_path = globals().get('analyses_path',os.path.join(base_path,'analyses'))
y_category = globals().get('y_category','Genotype')
boolean_types = globals().get('boolean_types',{'Genotype':{'KO':1,'WT':0},'Sex':{'Male':0,'Female':1}})
behavioural_tests = globals().get('behavioral_tests',["Novel Object Recognition (NOR) Training Preference","Novel Object Recognition (NOR) Testing Preference","NOR discrimination index"])
blood_tests = globals().get('blood_tests',["G-CSF/CSF-3","IL-10","IL-3","LIF","IL-1b","IL-2","M-CSF","IP-10/CXCL10","VEGF","IL-4","IL-5","IL-6","IL-25/IL-17","IFNa","IL-22","IL-9",
               "IL-13","IL-27","IL-23","IFNg","IL-12P70","GM-CSF","GROa/KC/CXCL1","RANTES/CCL5","TNF-a","sRANKL","MIP-1a/CCL3","MCP-3/CCL7","MCP-1/CCL2",
               "IL-17A/CTLA8","IL-7RA","IL-15/IL-15R","MIP-2","IL-1a","ENA-78/LIX/CXCL5","IL-19","EOTAXIN/CCL11","IL-2RA","IL-28","LEPTIN","IL-18",
               "BAFF/TNFSF13B","MIP-1b/CCL4","BTC","IL-7","IL-33","IL-31","ST2/IL-33R","CHEX-1","CHEX-2","CHEX-3","CHEX-4"])
time_behavioral = globals().get('time_behavioral','Age at behavioral tests (months)')
time_blood = globals().get('time_blood','Age at organ collection (months)')
time_dependence_plots =globals().get('time_dependence_plots',os.path.join(analyses_path,'time_dependence_plots'))