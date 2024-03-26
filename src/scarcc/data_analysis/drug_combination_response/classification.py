"""Classification of drug combination response based on predicted and observed values."""
# Current using po_col

import pandas as pd

def convert_po_col(col, additive_threshold=0.05):
    """Convert expected-observed column to classification based on threshold."""
    diff_bins = [-10,-1*additive_threshold,1*additive_threshold,10]
    return pd.cut(col, bins=diff_bins, labels=['Antagonistic', 'Additive', 'Synergistic'])

# TODO: replace po_col with op_col, flipping the sign and match Yeh paper observed-expected
def convert_op_col(col, additive_threshold=0.05):
    """Convert observed-expected column to classification based on threshold."""
    diff_bins = [-10,-1*additive_threshold,1*additive_threshold,10]
    return pd.cut(col, bins=diff_bins, labels=['Synergistic', 'Additive', 'Antagonistic'])
