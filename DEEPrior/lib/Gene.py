import sys
import os
import requests
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.Transcript import Transcript


class Gene:
    def __init__(self, obj_dict):
        self.start = obj_dict['start']
        self.end = obj_dict['end']
        self.ensg = obj_dict['gene_id']
        self.common_name = obj_dict['external_name']
        self.strand = obj_dict['strand']
        self.chr = obj_dict['seq_region_name']
        self.version = obj_dict['assembly_name']
        self.description = obj_dict['description']
        self.biotype = obj_dict['biotype']
        self.decoded = []
        self.transcripts = []

        self._transcript_query()
        if len(self.decoded) != 0:
            self._add_transcripts()

    def _transcript_query(self):
        """
        This method retrieves all TRANSCRIPTS, EXONS and CDS regions from ENSG name
        """
        if self.version == 'GRCh37':
            server = "https://grch37.rest.ensembl.org"
        elif self.version == 'GRCh38':
            server = "https://rest.ensembl.org"
        else:
            print('Please check ENSEMBL version, only grch37 or grch38 are allowed')
            sys.exit()

        ext = "/overlap/id/%s?feature=transcript;feature=exon;feature=CDS" % self.ensg

        r = requests.get(server + ext, headers={"Content-Type": "application/json"})

        if r.ok:
            self.decoded = r.json()

        return

    def _add_transcripts(self):
        # first select all available transcripts from self.decoded
        decod_transc = list(filter(lambda x: x['feature_type'] == 'transcript', self.decoded))

        transc_list = []
        for i in range(len(decod_transc)):
            elem = decod_transc[i]
            # check if parent ensg name is the desired one and create Transcript object only if it is protein coding
            # with CCDSID
            # If so, append the created Transcript object and append it to transc_list
            if elem['Parent'] == self.ensg and 'ccdsid' in elem.keys():
                if elem['biotype'] == 'protein_coding':
                    transc_list.append(Transcript(elem))

        decod_exon = list(filter(lambda x: x['feature_type'] == 'exon', self.decoded))

        ensts = list(map(lambda x: x.enst, transc_list))
        for i in decod_exon:
            # select only exons of transcripts identified with the previous loop and add to the corresponding transcript
            # object
            if i['Parent'] in ensts:
                ind = ensts.index(i['Parent'])
                transc_list[ind].add_exon(i)

        decod_cds = list(filter(lambda x: x['feature_type'] == 'cds', self.decoded))
        for i in decod_cds:
            # select only cds of transcripts identified with the previous loop and add to the corresponding transcript
            # object
            if i['Parent'] in ensts:
                ind = ensts.index(i['Parent'])
                transc_list[ind].add_cds(i)

        self.transcripts = transc_list

        return
