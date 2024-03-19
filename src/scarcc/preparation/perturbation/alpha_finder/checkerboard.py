import numpy as np
import pandas as pd
import itertools
import os

from .monoculture import MonocultureAlphaFinder

class CheckerboardAlphaFinder():
    def __init__(self, response_record, data_directory = None, **maf_kwargs) -> None:
        self.response_record = response_record
        self.target_gene_list = list(response_record.keys())
        self.species_list = list(response_record[self.target_gene_list[0]].keys())
        self.n_levels = len(self.response_record[self.target_gene_list[0]][self.species_list[0]]['nbiomass_x'])
        self.checker_lvs = self.construct_checker_lvs()
        self.result_dict = {}
        self.maf_kwargs = maf_kwargs
        self.data_directory = data_directory
        if 'acceptance_threshold_upper' not in self.maf_kwargs:
            self.maf_kwargs['acceptance_threshold_upper'] = .995
        if 'acceptance_threshold_lower' not in self.maf_kwargs:
            self.maf_kwargs['acceptance_threshold_lower'] = 1.001
        if 'precision' not in self.maf_kwargs:
            self.maf_kwargs['precision'] = 3
        print(self.n_levels)
    @staticmethod
    def get_new_response_record(**kwargs):
        maf = MonocultureAlphaFinder(**kwargs)
        maf.find_feasible_alpha()
        return maf

    def fill_dict(self, model, current_gene, nbiomass_x):
        nbiomass_lv = float(format(nbiomass_x, '.2f'))
        inner_dict = self.response_record[current_gene][model].setdefault('response', {})
        visited= inner_dict.keys()
        if float(format(nbiomass_lv, '.2f')) not in visited:
            if inner_dict.get(nbiomass_lv):
                closest_val = min(visited, key=lambda x: abs(x - nbiomass_x))
                search_alpha = self.response_record[current_gene][model]['response'][closest_val]
                search_alpha = search_alpha*1.1 if closest_val > nbiomass_x else search_alpha*0.9

            else:
                search_alpha = 1.2
            maf = self.get_new_response_record(
                model=model, search_alpha=search_alpha,
                current_gene=current_gene, target_normalized_biomass=nbiomass_lv,
                response_record=self.response_record, **self.maf_kwargs)
            return maf
        return self.response_record
    
    def construct_checker_lvs(self):
        checker_lvs = pd.DataFrame()
        MIC_levels = np.arange(0, self.n_levels)

        checker_lvs['lv_pairs'] = list(itertools.product(MIC_levels, MIC_levels))
        checker_lvs['lv1'] = checker_lvs['lv_pairs'].apply(lambda x: x[0])
        checker_lvs['lv2'] = checker_lvs['lv_pairs'].apply(lambda x: x[1])
        return checker_lvs
    
    def match_alpha_lvs(self, current_gene):
        species_alpha_table = {}
        for Species in self.species_list:
            record_in_gene_species = self.response_record[current_gene][Species]
            result_dict = {
            ith: {
                f'alpha': record_in_gene_species['response'][float(format(nbiomass, '.2f'))]['search_alpha'],
                'lv': ith,
                'Gene_inhibition': current_gene,
                'normalized_biomass': nbiomass}
                    for ith, nbiomass in enumerate(record_in_gene_species['nbiomass_x'])}
            species_alpha_table[Species] = pd.DataFrame.from_dict(result_dict, orient='index').set_index('Gene_inhibition').add_prefix(f'{Species.id}_')
        return species_alpha_table
    
    def process_record_in_gene(self, current_gene):
        outer_dict = self.match_alpha_lvs(current_gene)
        for sp1, sp2 in itertools.combinations(outer_dict.keys(), 2):
            self.result_dict[current_gene, sp1, sp2] = (
                self.checker_lvs.merge(outer_dict[sp1], left_on='lv1', right_on=f'{sp1.id}_lv')
                .merge(outer_dict[sp2], left_on=f'lv2', right_on=f'{sp2.id}_lv', suffixes=(f'_{sp1}', f'_{sp2}'))
                .drop(columns=['lv1', 'lv2']))
            self.result_dict[current_gene, sp1, sp2]['Gene_inhibition'] = current_gene
            self.result_dict[current_gene, sp1, sp2].set_index('Gene_inhibition', inplace=True)

    def get_checkerboard_alpha_table(self):
        for current_gene, sub_dict in self.response_record.items():
            for model, ss_dict in sub_dict.items():
                for nbiomass in sorted(ss_dict['nbiomass_x']):
                    AF = self.fill_dict(model, current_gene, nbiomass)
                    self.response_record = AF.response_record  if isinstance(AF, MonocultureAlphaFinder) else AF
                    alpha_lvs = ss_dict.setdefault('alpha_lvs', {})
                    IC_lv = float(format(nbiomass, '.2f'))
                    alpha_lvs.update({nbiomass: ss_dict['response'][IC_lv]['search_alpha']})

        [self.process_record_in_gene(current_gene) for current_gene in self.target_gene_list]
        self.alpha_table = pd.concat(self.result_dict.values(), axis=0).sort_values('lv_pairs')

        if self.data_directory is not None:
            self.alpha_table.to_csv(os.path.join(self.data_directory, 'alpha_table_checkerboard.csv'))
            print(f'Checkerboard alpha table saved to {self.data_directory}.')
        return self.alpha_table