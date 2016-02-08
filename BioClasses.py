class Gene(object):
    def __init__(self, source, pFrom, pTo, gid=None, organism=None, strand=None, cogid=None, arcogid=None):
        self.src = source
        self.gid = gid
        self.pFrom, self.pTo = pFrom, pTo
        self.pFrom = int(self.pFrom)
        self.pTo = int(self.pTo)
        self.strand = strand
        self.organism = organism
        self.cogid = cogid
        self.arcogid = arcogid
        self.tag=None

    def __str__(self):
        return  "%s_%d_%d"%(self.gid, self.pFrom,self.pTo)

    def set_sequence(self, gnmSeq, flank=0):
        self.seq = gnmSeq[self.pFrom-1-flank:self.pTo-1+flank]

    def get_translate(self):
        if self.strand=='+':
            return self.seq.translate()
        else:
            return self.seq.reverse_complement().translate()

    def overlaps(self, other):
        dist1 = other.pFrom - self.pTo
        dist2 = self.pFrom - other.pTo
        return True if dist1*dist2>=0 else False

    def overlap_distance(self, other):

        if self.overlaps(other):
            dist1 = self.pTo - other.pFrom
            dist2 = other.pTo - self.pFrom

            return dist1 if dist1>0 else dist2
        else:
            return 0

    def in_vicinity(self, other, threshold=3000):
        # if used massively, change it to avoid getting the distances twice,
        # first in overlaps and then in in_vicinity

        if self.overlaps(other):
            return True
        dist = other.pFrom - self.pTo if (other.pFrom - self.pTo)>0 else self.pFrom - other.pTo
        return True if dist<=threshold else False

    def set_cog(self, cogid):
        self.cogid = cogid

    def set_arcogid(self, arcogid):
        self.arcogid = arcogid

    def __cmp__(self, other):
        if self.pFrom>other.pFrom:
            return 1
        elif self.pFrom<other.pFrom:
            return -1
        else:
            return 0

    def set_product_enzyme(self, enzyme_type_name):
        self.enzyme_type_name = enzyme_type_name


class CrisprGene(object):
    def __init__(self, source, pFrom, pTo, gid=None, organism=None, strand=None, cogid=None,
                 cas_status=None, gene_name=None, is_seed=None, crispr_type=None, is_complete=None):

        self.src = source
        self.gid = gid
        self.pFrom, self.pTo = pFrom, pTo
        self.pFrom = int(self.pFrom)
        self.pTo = int(self.pTo)
        self.strand = strand
        self.organism = organism
        self.cogid = cogid
        self.cas_status = cas_status
        self.gene_name = gene_name
        self.is_seed = is_seed
        self.crispr_type = crispr_type
        self.is_complete = is_complete


    def __str__(self):
        return  "%s_%d_%d"%(self.gid, self.pFrom,self.pTo)

    def __cmp__(self, other):
        if self.pFrom>other.pFrom:
            return 1
        elif self.pFrom<other.pFrom:
            return -1
        else:
            return 0
