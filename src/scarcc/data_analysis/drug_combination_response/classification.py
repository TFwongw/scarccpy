import pandas as pd

# ! op_col
def convert_po_col(col, additive_threshold=0.05):
    diff_bins = [-10,-1*additive_threshold,1*additive_threshold,10] 
    return pd.cut(col, bins=diff_bins, labels=['Antagonistic', 'Additive', 'Synergistic'])
