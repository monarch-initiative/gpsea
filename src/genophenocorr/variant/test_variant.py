from ._variant import VariantCoordinates, verify_start_end_coordinates


@pytest.mark.parameterize('contig, start, end, ref, alt, chlen, expected',
                          ['1', 100, 101, 'G', 'C', 0, '1:101:101/C'],  # SNP

                          ['1', 100, 102, 'GG', 'G', -1, '1:101:102/G'],  # DEL
                          ['16', 89_284_128, 89_284_134, 'CTTTTT', 'C', 5, '16:89284129-89284134/C'],  # DEL
                          ['16', 89_284_087, 89_284_089, 'AC', 'A', -1, '16:89284088-89284089/A'],  # DEL

                          ['1', 100, 101, 'G', 'GC', 1, '1:102-101/C'],  # INS
                          ['16', 89_283_999, 89_284_000, 'A', 'AT', 1, '16:89284001:89284000/T'],  # custom INS

                          ['1', 100, 120, 'G', 'C', 0, 'TODO'],  # MNV
                          ['16', 89_283_999, 89_284_002, 'AAT', 'AGCG', 1, '16:89284000:89284002/AGCG'],  # custom MNV

                          # symbolic DEL
                          ['1', 100, 120, 'N', '<DEL>', -20, '1:101-120/DEL'],
                          ['9', 133_359_999, 133_360_011, 'N', '<DEL>', -12, '9:133360000-133360011/DEL'],

                          # symbolic DUP
                          ['1', 100, 103, 'N', '<DUP>', 3, '1:101-103/DUP'],
                          ['9', 133_359_999, 133_360_002, 'N', '<DUP>', 3, '9:133360000-133360002/DUP'],

                          # symbolic INS
                          # We can describe an insertion of 20 bases in two ways, depending on if
                          # the variant coordinates are trimmed. In the untrimmed case,
                          # the INS replaces the REF allele, in the trimmed case it does not.
                          # Therefore, the change length of the untrimmed is 19 - one base (`N`)
                          # is gone and 20 bases are inserted. The change length of the trimmed is 20 - zero bases
                          # are gone and 20 bases are inserted.
                          # In practical terms, the `cdna_start` and `cdna_end` fields of the VEP response will change
                          # depending on the trimming status.
                          ['1', 100, 101, 'N', '<INS>', 19, '1:102:101/INS'],
                          ['1', 101, 101, '',  '<INS>', 20, '1:102:101/INS'],

                          ['9', 133_359_999, 133_360_000, 'N', '<INS>', 29, '9:133360001-133360000/INS'],
                          ['9', 133_360_000, 133_360_000, '', '<INS>', 30, '9:133360001-133360000/INS'],
                          )
def test_verify_start_end_coordinates(contig, start, end, ref, alt, chlen, expected):
    vc = VariantCoordinates(contig, start, end, ref, alt, chlen)
    out = verify_start_end_coordinates(vc)
    assert out == expected
