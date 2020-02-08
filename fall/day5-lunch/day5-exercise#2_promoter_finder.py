#!/usr/bin/env Python3

"""
Usage: ./promoter_finder.py <ctab>
"""

import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep='\t', index_col='chr')

relevant_data = df.loc[:,['t_name','start', 'end','strand']]

relevant_data['promoter_start'] = relevant_data.loc[:,'start']
relevant_data['promoter_end'] = relevant_data.loc[:,'start']

forward_strand_filter = relevant_data.loc[:,'strand'] == '+'
relevant_data.loc[forward_strand_filter,'promoter_start'] = relevant_data.loc[forward_strand_filter,'start'] - 500
relevant_data.loc[forward_strand_filter,'promoter_end'] = relevant_data.loc[forward_strand_filter,'start'] + 500


reverse_strand_filter = relevant_data.loc[:,'strand'] == '-'
relevant_data.loc[reverse_strand_filter,'promoter_start'] = relevant_data.loc[reverse_strand_filter,'end'] - 500
relevant_data.loc[reverse_strand_filter,'promoter_end'] = relevant_data.loc[reverse_strand_filter,'end'] + 500


relevant_data = relevant_data.drop(columns=['start','end','strand'])
relevant_data = relevant_data[['promoter_start','promoter_end','t_name']]

filter_df = relevant_data.loc[:,'promoter_start'] < 1
relevant_data.loc[filter_df,'promoter_start'] = 1
relevant_data = relevant_data[['promoter_start','promoter_end','t_name']]


relevant_data.to_csv(sys.argv[2], sep='\t', header=False)