import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def get_row_color_legend():
    color_label = dict()
    color_label['Drug Combination Effect'] = {'Antagonistic': '#90fda9', 'Additive': '#ffff84', 'Synergistic': '#fed0fc', '': 'white'} # place holder
    color_label['rel_abund'] = {'E Dominance (0.7 ≤ Rel. Abund. < 0.75)': '#4622fc', 'E Dominance (0.5 ≤ Rel. Abund. < 0.7)': '#9d22fc', 
                                'S Dominance (Rel. Abund. > 0.5 for S)': '#da52d5', 'No Growth': 'grey'}
    color_label['all_label'] = {}
    color_label['all_label'].update(color_label['Drug Combination Effect'])
    color_label['all_label'].update(color_label['rel_abund'])

    patches = [mpatches.Patch(color=color, label=label) for  label, color in color_label['all_label'].items()]

    plt.figure(figsize=(8,1.46))
    plt.legend(handles=patches, loc='upper left', title='Drug Combination Effect       Species Relative Abundance (Rel. Abund.)      ', ncol=2,frameon=False)
    plt.gca().set_axis_off()
    plt.show()
