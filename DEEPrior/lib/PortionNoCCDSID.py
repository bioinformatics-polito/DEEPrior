import sys
import os
import requests
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.GeneNoCCDSID import GeneNoCCDSID


class PortionNoCCDSID:
    def __init__(self, chr_, coord, p, version='grch38'):
        self.p = p
        self.version = version
        self.breakpoint_chr = chr_
        self.breakpoint_coord = int(coord)
        self.breakpoint_info = []
        self.sequences = []
        self.genes = []
        self.ensg = ''
        self.enst = []
        self.genomic_regions = []
        self.common_name = ''
        self.biotype = ''

        self._gene_query()  # fills genes list with all info about the genes: ccdsid transcript(s), cds and exons

        # only if gene query has produced a valid output
        if len(self.genes) != 0:
            # only if there is at least one ccdsid transcript
            if len(self.genes[0].transcripts) != 0:
                # find breakpoint info
                self._breakpoint_info()
                # decide genomic_region for each transcript
                self._genomic_region()
                # retrieve sequence for each transcript
                self._retrieve_sequences()

    def _gene_query(self):
        """""
        This function performs an API query to ENSEMBL and returns the decoded list with all genes located in a
        specific GENOMIC region
        :param self, see __init__ function
        :return: self.genes
        """
        if self.version == 'grch37':
            server = "https://%s.rest.ensembl.org" % self.version
        elif self.version == 'grch38':
            server = "https://rest.ensembl.org"
        else:
            print('Please check ENSEMBL version, only grch37 or grch38 are allowed')
            sys.exit()

        feature = 'feature=gene'
        ext = "/overlap/region/human/%s:%s-%s?%s" % (self.breakpoint_chr, self.breakpoint_coord,
                                                     self.breakpoint_coord, feature)

        r = requests.get(server + ext, headers={"Content-Type": "application/json"})

        if r.ok:
            decoded = r.json()

            # only if decoded has a content
            if len(decoded) != 0:
                self.genes.append(GeneNoCCDSID(decoded[0]))
                self.ensg = self.genes[0].ensg
                self.common_name = self.genes[0].common_name
                self.biotype = self.genes[0].biotype

                # only if coding ccdsid transcripts are present
                if len(self.genes[0].transcripts) != 0:
                    for j in self.genes[0].transcripts:
                        self.enst.append(j.enst)
        return

    def _breakpoint_info(self):
        """
        This function returns in self.breakpoint_info the position of the breakpoint: UTR, CDS or INTRON
        Breakpoint info is returned for each transcript
        :return:
        """
        transcripts = self.genes[0].transcripts
        # for each transcript
        for i in range(len(transcripts)):
            exons = transcripts[i].exons_pos
            cds = transcripts[i].cds_pos
            # defines ranges of 5p and 3p UTR regions, depending on the strand
            # utr_ranges = [min(exons), min(cds), max(cds), max(exons)]
            if transcripts[i].strand == 1:
                utr_ranges_5p = [self.genes[0].start, min(cds)]
                utr_ranges_3p = [max(cds), self.genes[0].end]

            else:
                utr_ranges_5p = [max(cds), self.genes[0].end]
                utr_ranges_3p = [self.genes[0].start, min(cds)]

            # check if breakpoint is 5p UTR or 3p UTR
            if utr_ranges_5p[0] <= self.breakpoint_coord <= utr_ranges_5p[1]:
                info = 'UTR_5p'

            elif utr_ranges_3p[0] <= self.breakpoint_coord <= utr_ranges_3p[1]:
                info = 'UTR_3p'

            # if (utr_ranges[0] <= self.breakpoint_coord < utr_ranges[1]) or \
            #         (utr_ranges[2] < self.breakpoint_coord <= utr_ranges[3]):
            #     info = 'UTR'
            # check if breakpoint is cds or intron
            else:
                is_in_cds_range = []
                for j in range(0, len(cds), 2):
                    if cds[j] <= self.breakpoint_coord <= cds[j+1]:
                        is_in_cds_range.append(1)
                    else:
                        is_in_cds_range.append(0)

                if sum(is_in_cds_range) == 0:
                    info = 'intron'
                else:
                    info = 'cds'

            self.breakpoint_info.append(info)

        return

    def _genomic_region(self):
        """
        This method defines the genomic boundaries of the portion involved in the fusion, depending if the region
        involved is 5p or 3p
        :return:
        """
        for i in self.genes[0].transcripts:
            # if 5 prime
            if self.p == 5:
                if self.genes[0].strand == 1:
                    genomic = [min(i.exons_pos), self.breakpoint_coord]
                else:
                    genomic = [self.breakpoint_coord, max(i.exons_pos)]
            # if 3 prime
            elif self.p == 3:
                if self.genes[0].strand == 1:
                    genomic = [self.breakpoint_coord, max(i.exons_pos)]
                else:
                    genomic = [min(i.exons_pos), self.breakpoint_coord]

            self.genomic_regions.append(genomic)
        return

    def _retrieve_sequences(self):
        """
        This function is intended to build the final sequence
        :return: self.sequence
        """
        for i in range(len(self.genes[0].transcripts)):
            cds = self.genes[0].transcripts[i].cds_pos
            sorted_cds = sorted(cds)
            chr_ = self.breakpoint_chr
            query = []

            # determines the ordered genomic region of interest
            genomic_region = self.genomic_regions[i]
            sequence = ''

            # case beginning is in UTR, put UTR string in sequence
            if genomic_region[1] < sorted_cds[0]:
                sequence = 'UTR'

            # case the end is in UTR, put UTR string in sequence
            elif genomic_region[0] > sorted_cds[-1]:
                sequence = 'UTR'

            else:
                # set start and end if external to cds region
                if genomic_region[0] < sorted_cds[0]:
                    genomic_region[0] = sorted_cds[0]

                if genomic_region[1] > sorted_cds[-1]:
                    genomic_region[1] = sorted_cds[-1]

                add_region = False

                # for each piece of cds, check beginning and initialize add_region
                for k in range(0, len(sorted_cds)-1, 2):
                    # check if beginning is in cds
                    if sorted_cds[k] <= genomic_region[0] <= sorted_cds[k+1]:
                        stop = sorted_cds[k+1]
                        if genomic_region[1] < sorted_cds[k+1]:
                            stop = genomic_region[1]
                        query.append('"%s:%s..%s:1"' % (chr_, genomic_region[0], stop))

                        add_region = True

                    # only if there is another cds
                    if k < (len(sorted_cds) - 2):

                        # check if beginning is in intron
                        if sorted_cds[k+1] < genomic_region[0] < sorted_cds[k+2]:
                            stop = sorted_cds[k+2] - 1
                            if genomic_region[1] < sorted_cds[k + 2]:
                                stop = genomic_region[1]
                            query.append('"%s:%s..%s:1"' % (chr_, genomic_region[0], stop))

                            add_region = True

                        # add all cds if add_region has been previously set
                        if add_region:
                            # check the end
                            if genomic_region[1] >= sorted_cds[k+3]:
                                stop = sorted_cds[k+3]
                                query.append('"%s:%s..%s:1"' % (chr_, sorted_cds[k+2], stop))

                            elif sorted_cds[k+2] < genomic_region[1] < sorted_cds[k+3]:
                                stop = genomic_region[1]
                                query.append('"%s:%s..%s:1"' % (chr_, sorted_cds[k + 2], stop))

                # from query to sequence
                # request_post_ensembl allows a maximum of 50 post!!
                req = list(PortionNoCCDSID.chunks(query, 50))
                sequence = ''
                for j in req:
                    a = PortionNoCCDSID.request_post_ensembl(', '.join(j), 'grch37')
                    if a != '':
                        sequence += a
                    else:
                        sequence = ''
                        break

                # if the transcript is onto the reverse strand:
                if self.genes[0].strand == -1:
                    dna = Seq(sequence, generic_dna)
                    sequence = str(dna.reverse_complement())

            self.sequences.append(sequence)

        return

    @staticmethod
    def chunks(l, n):
        # Create a function called "chunks" with two arguments, l and n:
        # For item i in a range that is a length of l,
        for i in range(0, len(l), n):
            # Create an index range for l of n items:
            yield l[i:i+n]

    @staticmethod
    def request_post_ensembl(query_string, version='grch38'):
        """
        This method can accept a maximum of 50 post sequences
        :param query_string: MUST be in the format: '"X:1000000..1000100:1", "ABBA01004489.1:1..100"'
        :param version: grch37 or grch38. grch38 is the default parameter
        :return: dna sequence in the requested DNA regions
        """
        if version == 'grch37':
            server = "https://grch37.rest.ensembl.org"
        elif version == 'grch38':
            server = "https://rest.ensembl.org"
        else:
            print('Please check ENSEMBL version, only grch37 or grch38 are allowed')
            sys.exit()

        ext = "/sequence/region/human"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        r = requests.post(server + ext, headers=headers,
                          data='{ "regions" : [' + query_string + '] }')

        if r.ok:
            decoded = r.json()
            sequence = ''
            for i in decoded:
                sequence += i['seq']

        else:
            sequence = ''

        return sequence
