import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.Portion import Portion


class Fusion:
    def __init__(self, chr5p, coord5p, chr3p, coord3p, tissue, version='grch38'):
        self.chr5p = chr5p.replace('chr', '')
        self.chr3p = chr3p.replace('chr', '')
        self.coord5p = coord5p
        self.coord3p = coord3p
        self.tissue = tissue
        self.version = version
        self.portions = [Portion(self.chr5p, self.coord5p, 5, version), Portion(self.chr3p, self.coord3p, 3, version)]
        self.fusion_pair = ''
        self.info_breakpoints = []
        self.sequences = []
        self.protein_cod = []

        self._fusion_info()
        # only if both are protein coding, calculates info_breakpoints and sequences
        if all(elem == 'protein_coding' for elem in self.protein_cod):
            self._calculate_sequences()

    def _fusion_info(self):
        # calculate fusion_pair, protein_coding only if a there is a valid ENSG gene
        if len(self.portions[0].genes) != 0 and len(self.portions[1].genes) != 0:
            self.fusion_pair = self.portions[0].common_name + '_' + self.portions[1].common_name
            self.protein_cod.append(self.portions[0].biotype)
            self.protein_cod.append(self.portions[1].biotype)
        return

    def _calculate_sequences(self):
        """
        Calculates all possible combinations of portion 5p and portion 3p. Both genes must have at least a sequence
        and 'UTR' sequences are changed with ''
        :return: info_breakpoints and sequences
        """
        sequences5p = self.portions[0].sequences
        sequences3p = self.portions[1].sequences
        info5p = self.portions[0].breakpoint_info
        info3p = self.portions[1].breakpoint_info

        # permutations are performed only if both Portions have a sequence
        if '' not in sequences3p and '' not in sequences5p:
            for i in range(len(sequences5p)):
                if sequences5p[i] == 'UTR':
                    sequences5p[i] = ''

                for j in range(len(sequences3p)):
                    if sequences3p[j] == 'UTR':
                        sequences3p[i] = ''

                    self.info_breakpoints.append([info5p[i], info3p[j]])
                    self.sequences.append(sequences5p[i] + sequences3p[j])

        return

