import pandas
#from ALSElinav.UtilsALS import getOrganism,DimentionalReductionMergedDataSets
import os
import numpy as np

from eranb import plot_PCoA
from eranb.Extract import extract

figurePath='/Users/daphna/Downloads'
# times={'AkkKEGG':[0,80,90,110],'CoHousing':[60,80,100,120,140],
# 'WTvsSOD1':[63,71,85,91,125,142],'WTvsSOD10':[40,60,80,100]}
treatment=[(None,None)]
other_params={'organism':'human','diamond':False,'index':'index'}
#other_params = {'MostAbundant':True,'thresholdKegg':20 / ((1e6) * 0.8)}
#other_params = {'apply':np.mean,'thresholdKegg':20 / ((5e6) * 0.8),'diamond':True,'index':'UserID'}
#thresholdKegg = other_params['thresholdKegg']
thresholdKegg=0.0005
if not 'index' in other_params:
    other_params['index']='uniting_key'
for figure in ['Human']:#'AkkKEGG','WTvsSOD1','Human','WTvsSOD10'
    #for treatment in [('Akkermansia','PBS')]:#('SH','CH')]:#('Akkermansia','PBS')]:#
    for treatment in [(None,None)]:#('Akkermansia','PBS'),('Akkermansia','PBS')]:#(None,None)]:
        for genotype in [('Control','ALS')]:#('WT','SOD1')('Control','ALS'),('WT','WT'),('SOD1','SOD1')]:#('WT','SOD1')]:#('WT','WT')
            for method in ['PCoA']:#,'TSNE']:
                for source in ['mpa']:
                    '''[('all','%s_%s_%s_%s_%sAll'%(method,
                                            figure,"".join(treatment),"".join(str(genotype)),source)),
                                      ('D>0','%s_%s_%s_%s_%sExcludeDay0'%(method,
                                            figure,"".join(str(genotype)),"".join(treatment),source))]+'''
                    for timepoint,figname in [('all','PLEASE')]:
                            # [('all', '%s_%s_%s_%s_%sAll' % (method,
                            #                                 figure, "".join(treatment), "".join(str(genotype)),
                            #                                 source)),
                            #  ('D>0', '%s_%s_%s_%s_%sExcludeDay0' % (method,
                            #                                         figure, "".join(str(genotype)), "".join(treatment),
                            #                                         source))] :
                        # [('all', '%s_%s_%s_%s_%sAll' % (method,
                        #                                       figure,"".join(str(genotype)),
                        #                                     "".join(str(treatment)), source)),
                        #      ('D>0', '%s_%s_%s_%s_%sExcludeDay0' % (method,
                        #                                               figure, "".join(str(genotype)),
                        #                                             "".join(str(treatment)), source))]+\
                        # [(t,'%s_%s_%s_%s_%sDay%s'%(method,
                        #                    figure,genotype,"".join(str(treatment)),source,str(t))) \
                        #                     for t in times[figure]]:
                        SOD1micePBS = extract(\
                            source=source, figure=figure, genotype=genotype[0], \
                            treatment=treatment[1], timepoint=timepoint,**other_params)
                        SOD1miceAkkermansia = extract(\
                            source=source, figure=figure, genotype=genotype[1], \
                            treatment=treatment[0], timepoint=timepoint,**other_params)


                        figurename=os.path.join(figurePath,figname+'.pdf')
                        excelname=os.path.join(figurePath,figname+'.xlsx')
                        plot_PCoA.plot(SOD1micePBS, SOD1miceAkkermansia, legend=['%s %s '%(genotype[0],treatment[1]) + str(timepoint),
                                                                  '%s %s '%(genotype[1],treatment[0]) + str(timepoint)],
                                                        threshold=thresholdKegg, figurename=figurename,
                                                        method=method,excelName=excelname,
                                                                indexName=other_params['index'])