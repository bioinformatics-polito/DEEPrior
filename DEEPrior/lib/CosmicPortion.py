import os
import sys
import requests
import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.Transcript import Transcript


class CosmicPortion:
    def __init__(self, enst, cdna_coord, o, inserted_seq='', version='grch38'):
        self.enst = enst
        self.common_name = ''
        self.cdna_coord = cdna_coord
        self.o = o
        self.sequence = inserted_seq
        self.version = version
        self.transcript_obj = ''
        self.breakpoint_chr = []
        self.breakpoint_coord = []
        self.breakpoint_info = []
        self.protein_coding = ''

        # if not an inserted sequence
        if self.sequence == '':
            # convert cdna coordinates into genomic coordinates (e.g. [1, 535] becomes chromosomes [10, 10] and coord
            # [6166414, 61665880])
            self._convert_cdna_coordinates()

            # this part only if info about breakpoint are complete
            if '' not in self.breakpoint_chr:
                # retrieve transcript, exons positions and cds positions for the involved transcript.
                # All these info are put in self.transcript_obj
                self._exons_and_cds_from_enst()

                # only if the transcript is protein coding
                if self.protein_coding == 'protein_coding':
                    # retrieve info breakpoint (intron/utr/cds)
                    self._breakpoint_info()

                    # retrieve sequence
                    self._retrieve_sequence()

    @staticmethod
    def from_cdna_to_genomic(enst, cdna_pos, version='grch38'):
        if version == 'grch37':
            server = "https://grch37.rest.ensembl.org"
        elif version == 'grch38':
            server = "https://rest.ensembl.org"
        else:
            print('Please check ENSEMBL version, only grch37 or grch38 are allowed')
            sys.exit()

        ext = "/map/cdna/%s/%s..%s?" % (enst, cdna_pos, cdna_pos)

        r = requests.get(server + ext, headers={"Content-Type": "application/json"})

        if not r.ok:
            chrom = ''
            coord = ''
            strand = ''
        else:
            decoded = r.json()
            chrom = decoded['mappings'][0]['seq_region_name']
            coord = decoded['mappings'][0]['start']
            strand = decoded['mappings'][0]['strand']

        return chrom, coord, strand

    def _convert_cdna_coordinates(self):
        """
        This method parse COSMIC cdna format to obtain a genomic coordinates.
        (e.g. [1, 535] becomes chromosomes [10, 10] and coord [6166414, 61665880])
        :return: self.breakpoint_chr and breakpoint_coord
        """
        # check each element of ['1', '2568456-58']
        for i in self.cdna_coord:
            # keep track of upstream and downstream to convert coordinates
            upstream = 0
            downstream = 0
            if '+' in i:
                cdna_pos, downstream = i.split('+')
            elif '-' in i:
                cdna_pos, upstream = i.split('-')
            else:
                cdna_pos = int(i)

            chrom, coord, strand = CosmicPortion.from_cdna_to_genomic(self.enst, cdna_pos, 'grch37')

            if chrom != '':
                if downstream != 0:
                    if strand == 1:
                        coord += int(downstream)
                    else:
                        coord -= int(downstream)

                if upstream != 0:
                    if strand == 1:
                        coord -= int(upstream)
                    else:
                        coord += int(upstream)

            self.breakpoint_chr.append(chrom)
            self.breakpoint_coord.append(coord)
        return

    def _correctness(self):
        """
        This method check coordinates provided by Cosmic. E.g. from Cosmic definition ['1', '535'] 535-th cDNA base
        should be the end of an exon; ['1', '535+5688'] means that the breakpoint is 5688 bases after 535 cdna base
        :return: value equal to 0 if correct, 1 otherwise
        """
        values = []
        exons = self.transcript_obj.exons_pos
        len_exons = [(exons[i] - exons[i - 1] + 1) for i in range(1, len(exons), 2)]
        for cdna in self.portion:
            cdna_ = int(re.split('\\+|-', cdna)[0])

            cum = 0
            value = 1
            if cdna_ == 1:
                value = 0
            else:
                for i in len_exons:
                    if cdna_ == cum + i:
                        value = 0
            values.append(value)

        self.is_correct = values
        return

    def _exons_and_cds_from_enst(self):
        """
        This method retrieves all EXONS and CDS regions from ENST transcript name creating a Transcript object
        """
        if self.version == 'grch37':
            server = "https://grch37.rest.ensembl.org"
        elif self.version == 'grch38':
            server = "https://rest.ensembl.org"
        else:
            print('Please check ENSEMBL version, only grch37 or grch38 are allowed')
            sys.exit()

        ext = "/overlap/id/%s?feature=transcript;feature=exon;feature=CDS" % self.enst

        r = requests.get(server + ext, headers={"Content-Type": "application/json"})

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()

        # create transcript object and common name
        trasc = list(filter(lambda x: x['feature_type'] == 'transcript' and x['transcript_id'] == self.enst, decoded))
        self.transcript_obj = Transcript(trasc[0])
        self.common_name = trasc[0]['external_name'].split('-')[0]
        self.protein_coding = trasc[0]['biotype']

        # selecting only exons of self.enst and adding to self.transcript_obj
        exons = list(filter(lambda x: x['feature_type'] == 'exon' and x['Parent'] == self.enst, decoded))
        for i in exons:
            self.transcript_obj.add_exon(i)

        # selecting only cds of self.enst and adding to self.transcript_obj
        cds = list(filter(lambda x: x['feature_type'] == 'cds' and x['Parent'] == self.enst, decoded))
        for i in cds:
            self.transcript_obj.add_cds(i)

        return

    def _breakpoint_info(self):
        """
        This function returns in self.breakpoint_info the position of the breakpoint: UTR, CDS or INTRON
        :return:
        """
        exons = self.transcript_obj.exons_pos
        cds = self.transcript_obj.cds_pos
        # defines ranges of UTR regions
        utr_ranges = [min(exons), min(cds), max(cds), max(exons)]

        for i in self.breakpoint_coord:
            # check if breakpoint is UTR
            if (utr_ranges[0] <= i <= utr_ranges[1]) or (utr_ranges[2] <= i <= utr_ranges[3]):
                info = 'UTR'
            # check if breakpoint is cds or intron
            else:
                is_in_cds_range = []
                for j in range(0, len(cds), 2):
                    if cds[j] <= i <= cds[j+1]:
                        is_in_cds_range.append(1)
                    else:
                        is_in_cds_range.append(0)

                if sum(is_in_cds_range) == 0:
                    info = 'intron'
                else:
                    info = 'cds'

            self.breakpoint_info.append(info)

        return

    def _retrieve_sequence(self):
        """
        This function is intended to build the final sequence
        :return: self.sequence
        """
        cds = self.transcript_obj.cds_pos
        sorted_cds = sorted(cds)
        chr_ = self.breakpoint_chr[0]
        query = []

        # determines the ordered genomic region of interest
        genomic_region = [min(self.breakpoint_coord), max(self.breakpoint_coord)]

        # case beginning is in UTR, put UTR string in sequence
        if genomic_region[1] < sorted_cds[0]:
            self.sequence = 'UTR'

        # case the end is in UTR, put UTR string in sequence
        elif genomic_region[0] > sorted_cds[-1]:
            self.sequence = 'UTR'

        else:
            # set start and end if external to cds region
            if genomic_region[0] < sorted_cds[0]:
                genomic_region[0] = sorted_cds[0]

            if genomic_region[1] > sorted_cds[-1]:
                genomic_region[1] = sorted_cds[-1]

            add_region = False

            # for each piece of cds, check beginning and initialize add_region
            for k in range(0, len(sorted_cds) - 1, 2):
                # check if beginning is in cds
                if sorted_cds[k] <= genomic_region[0] <= sorted_cds[k + 1]:
                    stop = sorted_cds[k + 1]
                    if genomic_region[1] < sorted_cds[k + 1]:
                        stop = genomic_region[1]
                    query.append('"%s:%s..%s:1"' % (chr_, genomic_region[0], stop))

                    add_region = True

                # only if there is another cds
                if k < (len(sorted_cds) - 2):

                    # check if beginning is in intron
                    if sorted_cds[k + 1] < genomic_region[0] < sorted_cds[k + 2]:
                        stop = sorted_cds[k + 2] - 1
                        if genomic_region[1] < sorted_cds[k + 2]:
                            stop = genomic_region[1]
                        query.append('"%s:%s..%s:1"' % (chr_, genomic_region[0], stop))

                        add_region = True

                    # add all cds if add_region has been previously set
                    if add_region:
                        # check the end
                        if genomic_region[1] >= sorted_cds[k + 3]:
                            stop = sorted_cds[k + 3]
                            query.append('"%s:%s..%s:1"' % (chr_, sorted_cds[k + 2], stop))

                        elif sorted_cds[k + 2] < genomic_region[1] < sorted_cds[k + 3]:
                            stop = genomic_region[1]
                            query.append('"%s:%s..%s:1"' % (chr_, sorted_cds[k + 2], stop))

            # from query to sequence
            # request_post_ensembl allows a maximum of 50 post!!
            req = list(CosmicPortion.chunks(query, 50))
            for i in req:
                self.sequence += CosmicPortion.request_post_ensembl(', '.join(i), 'grch37')

            # if the transcript is onto the reverse strand:
            if self.transcript_obj.strand == -1:
                dna = Seq(self.sequence, generic_dna)
                self.sequence = str(dna.reverse_complement())

            # if COSMIC report antisense 'o' orientation:
            if self.o == 1:
                antisense = Seq(self.sequence, generic_dna)
                self.sequence = str(antisense.reverse_complement())
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

        if not r.ok:
            r.raise_for_status()
            sys.exit()

        decoded = r.json()
        sequence = ''
        for i in decoded:
            sequence += i['seq']

        return sequence
