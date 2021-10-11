import sys
import os
import requests
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)))
from lib.Transcript import Transcript


class GeneNoCCDSID:
    def __init__(self, obj_dict):
        try:
            self.start = obj_dict['start']
        except KeyError:
            self.start = []
        try:
            self.end = obj_dict['end']
        except KeyError:
            self.end = []
        try:
            self.ensg = obj_dict['gene_id']
        except KeyError:
            self.ensg = []
        try:
            self.common_name = obj_dict['external_name']
        except KeyError:
            self.common_name = []
        try:
            self.strand = obj_dict['strand']
        except KeyError:
            self.strand = []
        try:
            self.chr = obj_dict['seq_region_name']
        except KeyError:
            self.chr = []
        try:
            self.version = obj_dict['assembly_name']
        except KeyError:
            self.version = []
        try:
            self.description = obj_dict['description']
        except KeyError:
            self.description = []
        try:
            self.biotype = obj_dict['biotype']
        except KeyError:
            self.biotype = []
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
            # check if parent ensg name is the desired one
            # Then create Transcript object and append it to transc_list
            if elem['Parent'] == self.ensg:
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
