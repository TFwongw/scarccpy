from scarcc.preparation.target_gene.gene_format_handler import GeneFormatHandler

def test_GeneFormatHandler():
    assert GeneFormatHandler('folP').DG is None
    assert GeneFormatHandler('folP').SG == 'folP'
    assert GeneFormatHandler('folP').is_checkerboard_format is None
    assert GeneFormatHandler('folP.folA').DG is not None
    assert GeneFormatHandler('folP.folA').first_gene == 'folP'
    assert GeneFormatHandler('folP.folA').second_gene == 'folA'
    assert GeneFormatHandler('folP.folA_0.1').DG is None
    assert GeneFormatHandler('folP.folA_0.1').is_checkerboard_format is True
    assert GeneFormatHandler('folP.folA_0.1').first_gene == 'folP_0'
    assert GeneFormatHandler('folP.folA_0.1').second_gene == 'folA_1'
    assert GeneFormatHandler('folP.folA_0.0').DG is None
    assert GeneFormatHandler('folP.folA_2.3').DG is not None
    assert GeneFormatHandler('folP.folA_2.3').first_gene == 'folP.folA_2.0'
    assert GeneFormatHandler('folP.folA_2.3').second_gene == 'folP.folA_0.3'

    assert GeneFormatHandler('S0_folP.folA_coculture_2.3').DG == 'folP.folA_2.3'
    assert GeneFormatHandler('E0_folP.folA_monoculture_2.3').second_gene == 'folP.folA_0.3'
    assert GeneFormatHandler('folP').SG == 'folP'
    assert GeneFormatHandler('E0_folP').SG == 'folP'

    assert GeneFormatHandler('E0_Normal').is_Normal is True
    assert GeneFormatHandler('E0_folP.folA_0.0').is_Normal is True
    assert GeneFormatHandler('E0_folP.folA_1.0').is_Normal is False
    assert GeneFormatHandler('E0_folP.folA_3.2').is_Normal is False