import re
import requests
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.CosmicPortion import CosmicPortion as Portion


class Cosmic:
    def __init__(self, translocation_name, tissue, version='grch38'):
        self.fusion_pair = ''
        self.fusion_enst = ''
        self.tissue = tissue
        self.portions = []  # list of Portion objects, containing the pieces of the fusion
        self.translocation_name = translocation_name
        self.version = version
        self.breakpoints = []
        self.chr_breakpoints = []
        self.info_breakpoints = []
        self.sequence = ''
        self.protein_cod = []

        self._fusion_info()  # fills self.portions, fusion_pair and fusion_enst

        # only if len(self.portions) >= 2 and portions.protein_coding are all coding

        if len(self.portions) >= 2 and '' not in [j.sequence for j in self.portions]:
            self._fusion_break()  # obtain 'standard' info about the fusion
            self._final_seq()  # fills self.sequence

    def _fusion_info(self):
        """
        This method is intended to parse COSMIC translocation names and to obtain the final sequence
        :return:  Cosmic fusion information
        """

        to_split = self.translocation_name.replace('NM_', 'NM')
        split = re.split('{|}:r.|_', to_split)

        # for each piece of the fusion:
        k = 0  # offset if 'ins' is present in translocation name
        for i in range(0, len(split)-2, 4):
            # check if piece is 'ins_' and, if true, check for ins validity
            if 'ins' in split[i+k]:
                valid_ins = Cosmic.ins_valid(split[i+k])
                # print(valid_ins)
                k += 1
                if valid_ins != -1:
                    self.portions.append(Portion('', '', 0, valid_ins, self.version))
                else:
                    self.portions = []
                    break

            # define enst, portion, reverse
            enst = split[i+k+1]
            if 'NM' in enst:
                enst = self._nm_query((enst.replace('NM', 'NM_')).split('.')[0])

            # check for reference sequences different than 'ENST....'. In case, break portion loop
            if 'ENST' not in enst:
                self.portions = []
                break

            portion = split[i+k+2:i+k+4]
            # check if 'r.' is in portion. In case break portion loop
            if 'r.' in portion[0] or 'r.' in portion[1]:
                self.portions = []
                break

            # check if reverse
            if split[i+k][0] == 'o':
                o = 1
            else:
                o = 0

            # add to self.portions a portion object
            self.portions.append(Portion(enst, portion, o, '', self.version))

        # add self.fusion_pair, self.fusion_enst, self.sequence, self.fusion_breakpoint
        if len(self.portions) >= 2:
            self.fusion_pair = self.portions[0].common_name + '_' + self.portions[-1].common_name
            self.fusion_enst = self.portions[0].enst + '_' + self.portions[-1].enst
            for i in self.portions:
                self.protein_cod.append(i.protein_coding)

        return

    def _fusion_break(self):
        """
        This method is intended to obtain standard info about a Cosmic Fusion (i.e. chromosome, coordinates, sequences)
        :return:
        """

        # first element of the fusion
        self.breakpoints.append(self.portions[0].breakpoint_coord[1])
        self.chr_breakpoints.append(self.portions[0].breakpoint_chr[1])
        self.info_breakpoints.append(self.portions[0].breakpoint_info[1])

        # last element of the fusion
        self.breakpoints.append(self.portions[-1].breakpoint_coord[0])
        self.chr_breakpoints.append(self.portions[-1].breakpoint_chr[0])
        self.info_breakpoints.append(self.portions[-1].breakpoint_info[0])
        return

    def _final_seq(self):
        """
        This function retrieves the overall sequence of a Cosmic object, concatenating all portion sequences
        :return:
        """
        for i in range(len(self.portions)):
            if self.portions[i].sequence == 'UTR':
                self.sequence += ''
            else:
                self.sequence += self.portions[i].sequence

        return

    def _nm_query(self, nm):
        """""
        This function performs an API query to ENSEMBL and returns the enst name of an NM_ gene
        :param self, see __init__ function
        :return: decoded list
        """
        if self.version == 'grch37':
            server = "https://%s.rest.ensembl.org" % self.version
        elif self.version == 'grch38':
            server = "https://rest.ensembl.org"
        else:
            print('Please check ENSEMBL version, only grch37 or grch38 are allowed')
            sys.exit()

        ext = "/xrefs/symbol/homo_sapiens/" + nm + "?"
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()
        enst = decoded[1]['id']

        return enst

    @staticmethod
    def ins_valid(ins):
        """
        This method check is inserted sequence is fine (e.g. 'insAUCGUAGC') or not (eg. 'ins34' or 'ins?')
        :param ins: the ins string to verify (e.g. 'insAUCGUAGC')
        :return: the correct string (e.g. 'ATCGTAGC') or -1 if ins sequence is not valid
        """
        seq = ins.replace('ins', '').replace('U', 'T')
        chars = ['A', 'T', 'C', 'G']
        if any((c not in chars) for c in seq):
            seq = -1
        return seq
