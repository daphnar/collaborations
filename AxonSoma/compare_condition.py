import pandas as pd
import os
os.environ['OAK_BASE']='/oak/stanford/groups/pritch'
output_base_folder = os.path.join(os.environ['OAK_BASE'],'users/daphna/rdna/analyses/nanopore/pipeline/AxonSoma28s/')
nuc_path = '/home/users/daphna/analyses/nanopore/pipeline/AxonSoma28s/%s/haplotypes/mouse_28s/rRNA.mouse_28s_chunk_0.haplotypes.nuc.csv'
interesting_threshold_change = 0.01
samples = ['m64283e_240201_004140.SomaGFP--SomaGFP','m64283e_240201_004140.AxonGFP--AxonGFP',
           'm64283e_240201_004140.AxonTDP--AxonTDP','m64283e_240201_004140.SomaTDP--SomaTDP',\
           'm64283e_240201_004140.SomaPR--SomaPR','m64283e_240201_004140.SomaNone--SomaNone']
variants_with_at_least_5_percent=[]
for sample in samples:
    cond = pd.read_csv(nuc_path % sample, index_col=0)
    for idx in range(4729):
        if (cond[str(idx)].value_counts()[1]) >= ((cond.shape[0] - 1) * 0.05):
            variants_with_at_least_5_percent.append(str(idx))

pre_calculated_variants_with_at_least_5_percent = [  11,   12,   13,   14,   15,   16,   17,   18,   19,   20,   21,
         22,   23,   24,   25,   26,   27,   28,   29,   30,   31,   32,
         33,   34,   35,   36,   37,   38,   39,   40,   41,   42,   43,
         44,   45,   46,   47,   48,   49,   50,   51,   52,   53,   54,
         55,   56,   57,   58,   59,   60,   61,   62,   63,   64,   65,
         66,   67,   68,   69,   70,   71,   72,   73,   74,   75,   76,
         77,   78,   79,   80,   81,   82,   83,   84,   85,   86,   87,
         88,   89,   90,   91,   92,   93,   94,   95,   96,   97,   98,
         99,  100,  101,  102,  103,  104,  105,  106,  107,  108,  109,
        110,  111,  112,  113,  114,  115,  116,  117,  118,  119,  120,
        121,  122,  123,  124,  125,  126,  127,  128,  129,  130,  131,
        132,  133,  134,  135,  136,  137,  138,  139,  140,  141,  142,
        143,  187,  601,  651,  762,  882,  961, 1135, 1512, 1891, 1892,
       1893, 1894, 1895, 2363, 2502, 2822, 3102, 3147, 3201, 3715, 4181,
       4482, 4613, 4614, 4615, 4616, 4617, 4618, 4619, 4620, 4621, 4622,
       4623, 4624, 4625, 4626, 4627, 4628, 4629, 4630, 4631, 4632, 4633,
       4634, 4635, 4636, 4637, 4638, 4639, 4640, 4641, 4642, 4643, 4644,
       4645, 4646, 4647, 4648, 4649, 4650, 4651, 4652, 4653, 4654, 4655,
       4656, 4657, 4658, 4659, 4660, 4661, 4662, 4663, 4664, 4665, 4666,
       4667, 4668, 4669, 4670, 4671, 4672, 4673, 4674, 4675, 4676, 4677,
       4678, 4679, 4680, 4681, 4682, 4683, 4684, 4685, 4686, 4687, 4688,
       4689, 4690, 4691, 4692, 4693, 4694, 4695, 4696, 4697, 4698, 4699,
       4700, 4701, 4702, 4703, 4704, 4705, 4706, 4707, 4708, 4709, 4710,
       4711, 4712, 4713, 4714, 4715, 4716]
pre_calculated_variants=[187, 601, 651, 762, 882, 961, 1135, 1512, 1891, 1892,
1893, 1894, 1895, 2363, 2502, 2822, 3102, 3147, 3201, 3715, 4181,
4482]
pre_calculated_variants_leading=pd.DataFrame(index=[187, 601, 651, 762, 882, 961, 1135, 1512, 1891, 1892,
1893, 1894, 1895, 2363, 2502, 2822, 3102, 3147, 3201, 3715, 4181,
4482],columns=['leading_nuc','m64283e_240201_004140.SomaGFP--SomaGFP','m64283e_240201_004140.AxonGFP--AxonGFP',
           'm64283e_240201_004140.AxonTDP--AxonTDP','m64283e_240201_004140.SomaTDP--SomaTDP',\
           'm64283e_240201_004140.SomaPR--SomaPR','m64283e_240201_004140.SomaNone--SomaNone'])
pre_calculated_variants_alt=pd.DataFrame(index=[187, 601, 651, 762, 882, 961, 1135, 1512, 1891, 1892,
1893, 1894, 1895, 2363, 2502, 2822, 3102, 3147, 3201, 3715, 4181,
4482],columns=['alt_nuc','m64283e_240201_004140.SomaGFP--SomaGFP','m64283e_240201_004140.AxonGFP--AxonGFP',
           'm64283e_240201_004140.AxonTDP--AxonTDP','m64283e_240201_004140.SomaTDP--SomaTDP',\
           'm64283e_240201_004140.SomaPR--SomaPR','m64283e_240201_004140.SomaNone--SomaNone'])

# conditions=[['m64283e_240201_004140.SomaPR--SomaPR','m64283e_240201_004140.SomaGFP--SomaGFP'],
#             ['m64283e_240201_004140.SomaPR--SomaPR','m64283e_240201_004140.SomaNone--SomaNone'],
#             ['m64283e_240201_004140.SomaTDP--SomaTDP','m64283e_240201_004140.SomaGFP--SomaGFP'],
#             ['m64283e_240201_004140.SomaTDP--SomaTDP','m64283e_240201_004140.SomaNone--SomaNone'],
#             ['m64283e_240201_004140.AxonTDP--AxonTDP','m64283e_240201_004140.AxonGFP--AxonGFP'],
#             ['m64283e_240201_004140.AxonTDP--AxonTDP','m64283e_240201_004140.SomaTDP--SomaTDP'],
#             ['m64283e_240201_004140.AxonGFP--AxonGFP','m64283e_240201_004140.SomaGFP--SomaGFP']]
#interesting_variants = {}
# variants_with_at_least_5_percent = []
# for condition_pair in conditions:
#     first_cond = pd.read_csv(nuc_path%condition_pair[0],index_col=0)
#     second_cond = pd.read_csv(nuc_path%condition_pair[1],index_col=0)
#     for idx in range(4729):
#         if (first_cond[str(idx)].value_counts()[1])>=((first_cond.shape[0]-1)*0.05) or \
#             (second_cond[str(idx)].value_counts()[1]) >= ((second_cond.shape[0]-1)*0.05):
#             variants_with_at_least_5_percent.append(str(idx))
pre_calculated_variants_leading['leading_nuc']=None
pre_calculated_variants_alt['alt_nuc']=None

for sample in samples:
    cond = pd.read_csv(nuc_path % sample, index_col=0)
    for position in pre_calculated_variants:
        leading_nuc = cond[str(position)].value_counts().index[0]
        alt_nuc = cond[str(position)].value_counts().index[1]
        if pre_calculated_variants_leading.loc[position,'leading_nuc'] is None:
            pre_calculated_variants_leading.loc[position,'leading_nuc']=leading_nuc
        else:
            assert(pre_calculated_variants_leading.loc[position,'leading_nuc']==leading_nuc)
        if pre_calculated_variants_alt.loc[position, 'alt_nuc'] is None:
            pre_calculated_variants_alt.loc[position,'alt_nuc']=alt_nuc
            assert (pre_calculated_variants_alt.loc[position, 'alt_nuc'] == alt_nuc)
        pre_calculated_variants_leading.loc[position, sample]=cond[str(position)].value_counts()[0]/cond.shape[0]
        pre_calculated_variants_alt.loc[position, sample]=cond[str(position)].value_counts()[1]/cond.shape[0]

pre_calculated_variants_leading.to_csv('/home/users/daphna/analyses/nanopore/pipeline/AxonSoma28s/leading_variant_changes.csv')
pre_calculated_variants_alt.to_csv('/home/users/daphna/analyses/nanopore/pipeline/AxonSoma28s/alt_variant_changes.csv')


