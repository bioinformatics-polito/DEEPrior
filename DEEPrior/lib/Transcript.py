class Transcript:
    def __init__(self, obj_dict):
        self.enst = obj_dict['transcript_id']
        self.biotype = obj_dict['biotype']
        self.chr_ = obj_dict['seq_region_name']
        self.strand = obj_dict['strand']
        self.start = obj_dict['start']
        self.end = obj_dict['end']
        self.exons_pos = []
        self.cds_pos = []

    def add_exon(self, ex_dict):
        self.exons_pos.append(ex_dict['start'])
        self.exons_pos.append(ex_dict['end'])

    def add_cds(self, cds_dict):
        self.cds_pos.append(cds_dict['start'])
        self.cds_pos.append(cds_dict['end'])

