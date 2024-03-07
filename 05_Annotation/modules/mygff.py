"""
Created on May 09, 2017

Version 1.0

@author: sven.twardziok@posteo.de
"""

import csv, re, math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

class Feature(object):
    """Class for single features
    
    based on gff3 specification:
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    """
    
    def __lt__(self, other):
        """Defines behavior for the less-than operator, <
        
        :param other: other feature object to compare with
        :type other: object
        """
        
        if self.seqid<other.seqid or (self.seqid==other.seqid and self.start<other.start):
            return True
        elif self.seqid==other.seqid and self.start==other.start and self.end<other.end:
            if self.ftype in ["gene"]:
                return True
            elif self.ftype in ["mRNA"] and other.ftype not in ["gene"]:
                return True
            elif self.ftype in ["exon"] and other.ftype not in ["mRNA", "gene"]:
                return True
            elif self.ftype in ["three_prime_UTR", "five_prime_UTR", "CDS", "intron"] and other.ftype not in ["exon", "mRNA", "gene"]:
                return True
        elif self.seqid==other.seqid and self.start==other.start and self.end>=other.end:
            if self.ftype in ["gene"] and other.ftype not in ["gene"]:
                return True
            elif self.ftype in ["mRNA"] and other.ftype not in ["mRNA", "gene"]:
                return True
            elif self.ftype in ["exon"] and other.ftype not in ["exon", "mRNA", "gene"]:
                return True
        else:
            return False
    
    def __gt__(self, other):
        """Defines behavior for the greater-than operator, >
        
        :param other: other feature object to compare with
        :type other: object
        """
        
        if self.seqid>other.seqid or (self.seqid==other.seqid and self.start>other.start):
            return True
        elif self.seqid==other.seqid and self.start==other.start and self.end>other.end:
            if self.ftype in ["gene"] and other.ftype in ["gene"]:
                return True
            if self.ftype in ["mRNA"] and other.ftype in ["mRNA", "gene"]:
                return True
            if self.ftype in ["exon"] and other.ftype in ["mRNA", "gene", "exon"]:
                return True
            if self.ftype in ["three_prime_UTR", "five_prime_UTR", "CDS", "intron"]:
                return True
        elif self.seqid==other.seqid and self.start==other.start and self.end<=other.end:
            if self.ftype in ["mRNA"] and other.ftype in ["gene"]:
                return True
            if self.ftype in ["exon"] and other.ftype in ["mRNA", "gene"]:
                return True
            if self.ftype in ["three_prime_UTR", "five_prime_UTR", "CDS", "intron"] and other.ftype in ["exon", "mRNA", "gene"]:
                return True
        else:
            return False
    
    def __eq__(self, other):
        """Defines behavior for the equality operator, ==
        
        :param other: other feature object to compare with
        :type other: object
        """
        
        if self.seqid==other.seqid and self.start==other.start and self.end==other.end:
            # define equality for mRNAs
            if self.ftype=="mRNA" and other.ftype=="mRNA":
                # get all CDSs from both mRNAs
                cdss_self = []
                cdss_other = []
                for feature in self.features:
                    if feature.ftype=="CDS":
                        cdss_self.append(feature)
                for feature in other.features:
                    if feature.ftype=="CDS":
                        cdss_other.append(feature)
                cdss_self = sorted(cdss_self)
                cdss_other = sorted(cdss_other)
                # check if number of CDSs are equal
                if len(cdss_self) == len(cdss_other):
                    # check if all CDSs are equal and return False if one pair is unequal
                    for i in range(0, len(cdss_self)):
                        if cdss_self[i] != cdss_other[i]:
                            return False
                    return True
            elif self.ftype==other.ftype:
                return True
        else:
            return False
    
    def __hash__(self):
        return hash((self.seqid, self.start, self.end, self.ftype, self.identifier))
    
    def __init__(self, seqid, source, ftype, start, end, score, strand, phase):
        """Create feature object
        
        :param seqid: sequence identifier
        :type seqid: string
        :param source: name of source
        :type source: string
        :param ftype: feature type ("exon", "mRNA", "gene", "three_prime_UTR", "five_prime_UTR", "CDS", "intron")
        :type ftype: string
        :param start: start position
        :type start: int
        :param end: end position
        :type end: int 
        :param score: score value
        :type score: imt
        :param strand: strand inforamtion
        :type strand: string ("+", "-" or ".")
        :param phase: phase information
        :type phase: string
        """
        
        # standard fields from gff3 columns
        self.seqid = seqid
        self.source = source
        self.ftype = ftype
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        
        # attributes
        self.identifier = ""
        self.name = ""
        self.alias = ""
        self.notes = ""
        self.target = ""

        
        #links between features
        self.parent = None
        self.features = []
        
        # annotation stuff
        self.primary_confidence_class = ""
        self.secondary_condidence_class = ""
        
    
    def getLine(self):
        writeattributes = ""
        # required attributes
        if self.ftype=="gene":
            writeattributes = "ID=%s" % (self.identifier)
        elif self.ftype=="mRNA":
            if self.parent is None:
                print("error, no parent for: %s" % (self.identifier))
            else:
                writeattributes = "ID=%s;Parent=%s" % (self.identifier, self.parent.identifier)
        else:
            if self.parent is None:
                print("error, no parent for: %s %s %s %i %i" % (self.seqid, self.source, self.ftype, self.start, self.end))
            else:
                writeattributes = "Parent=%s" % (self.parent.identifier)
        #optional attributes
        if len(self.name)>0:
            writeattributes += ";Name=%s" % (self.name)
        if len(self.alias)>0:
            writeattributes += ";Alias=%s" % (self.alias)
        if len(self.target)>0:
            writeattributes += ";Target=%s" % (self.target)
        if len(self.notes)>0:
            writeattributes += ";Notes=%s" % (self.notes)
        if len(self.primary_confidence_class)>0:
            writeattributes += ";primary_confidence_class=%s" % (self.primary_confidence_class)
        if len(self.secondary_condidence_class)>0:
            writeattributes += ";secondary_confidence_class=%s" % (self.secondary_condidence_class)
        
        return [self.seqid, self.source, self.ftype, self.start, self.end, self.score, self.strand, self.phase, writeattributes]


class GeneAnnotation(object):
    """Read specific gff files and returns structured data for plant.annot"""
    
    def readGff3PlantAnnot(self, path):
        """General GFF3 file used in plant.annot pipeline
        
        :param path: path to gff file
        :type path: string
        
        0 seqname chrX       Chromosome, scaffold or contig name
        1 source  name       Name of source, e.g. database or software
        2 feature exon       "three_prime_UTR", "five_prime_UTR", "mRNA", "exon", "CDS", "gene", "intron"
        3 start   77696957   The leftmost coordinate of this record (where 1 is the leftmost possible coordinate)
        4 end     77712009   The rightmost coordinate of this record, inclusive.
        5 score   0.3221     Some score value
        6 strand  +          One of "+", "-", "."
        7 frame   .          Frame for feature (just used for CDS)
        8 attributes (GFF3)  ID=XXX;Parent=XXX (ID is only used for genes and mRNAs; Parent is not used for genes)
        """
        
        self.features = []
        self.genes = {}
        self.mrnas = {}
        self.seqids = {}
        genes2mrnas = []
        mrnas2features = []
        
        # create features
        with open(path, "r") as ingff3:
            reader = csv.reader(ingff3, delimiter="\t", quoting = csv.QUOTE_NONE)
            for line in reader:
                if len(line)==9:
                    seqid = line[0]
                    source = line[1]
                    ftype = line[2]
                    start = int(line[3])
                    end = int(line[4])
                    score = line[5]
                    strand = line[6]
                    phase = line[7]
                    feature = Feature(seqid, source, ftype, start, end, score, strand, phase)
                    attributesline = line[8]
                    attributes = {}
                    for entry in attributesline.split(";"):
                        matchAttribute = re.match(r"(.*)=(.*)", entry)
                        if matchAttribute:
                            attributes[matchAttribute.group(1)] = matchAttribute.group(2)
                    # add attributes to feature
                    if "ID" in attributes.keys():
                        feature.identifier = attributes["ID"]
                    if "Name" in attributes.keys():
                        feature.name = attributes["Name"]
                    if "Alias" in attributes.keys():
                        feature.alias = attributes["Alias"]
                    if "Notes" in attributes.keys():
                        feature.notes = attributes["Notes"]
                    if "Target" in attributes.keys():
                        feature.target = attributes["Target"]
                    if "primary_confidence_class" in attributes.keys():
                        feature.primary_confidence_class = attributes["primary_confidence_class"]
                    if "secondary_condidence_class" in attributes.keys():
                        feature.secondary_condidence_class = attributes["secondary_condidence_class"]
                    if "primconf" in attributes.keys():
                        feature.primary_confidence_class = attributes["primconf"] #old version
                    if "secconf" in attributes.keys():
                        feature.secondary_condidence_class = attributes["secconf"] #old version
                    # add gene to seqid and genes
                    if feature.ftype == "gene":
                        self.features.append(feature)
                        if not feature.seqid in self.seqids.keys():
                            self.seqids[seqid] = []
                        self.seqids[feature.seqid].append(feature)
                        self.genes[feature.identifier] = feature
                    # add mrna to mrnas and mark for gene assignment
                    elif feature.ftype == "mRNA":
                        self.features.append(feature)
                        self.mrnas[feature.identifier] = feature
                        genes2mrnas.append({"geneid":attributes["Parent"], "mrna":feature})
                    # mark remaining features for mrna assignment
                    elif feature.ftype in ["exon", "three_prime_UTR", "five_prime_UTR", "CDS", "intron"]:
                        self.features.append(feature)
                        mrnas2features.append({"mrnaid":attributes["Parent"], "feature":feature})
        
        # assign genes to mrnas
        for assignment in genes2mrnas:
            geneid = assignment["geneid"]
            mrna = assignment["mrna"]
            if geneid in self.genes.keys():
                gene = self.genes[geneid]
                mrna.parent = gene
                gene.features.append(mrna)
            else:
                print("gene missing")
        
        # assign mrnas to features
        for assignment in mrnas2features:
            mrnaid = assignment["mrnaid"]
            feature = assignment["feature"]
            if mrnaid in self.mrnas.keys():
                mrna = self.mrnas[mrnaid]
                feature.parent = mrna
                mrna.features.append(feature)
            else:
                print("mrna missing")
        
        # sort features and return
        self.features = sorted(self.features)
        return(self)
    
    def combine(self, geneannotations, annoversion="PGSB"):
        self.features = []
        for geneannotation in geneannotations:
            self.features += geneannotation.features
        self.features = sorted(self.features)
        genecounter = 0
        mrnacounter = 0
        self.genes = {}
        self.mrnas = {}
        self.seqids = {}
        for feature in self.features:
            if feature.ftype=="gene":
                genecounter += 1
                if not feature.seqid in self.seqids.keys():
                    self.seqids[feature.seqid] = []
                feature.identifier = "%s_gene_%i" % (annoversion, genecounter)
                self.genes[feature.identifier] = feature
                self.seqids[feature.seqid].append(feature)
            if feature.ftype=="mRNA":
                mrnacounter += 1
                feature.identifier = "%s_mRNA_%i" % (annoversion, mrnacounter)
                self.mrnas[feature.identifier] = feature
        return(self)
    
    def recalcGeneids(self, annoversion="PGSB"):
        #1) get one new gene for each mrna; import all attributes from former genes
        tmpnewgenes = []
        tmpcounter = 0
        for feature in self.features:
            if feature.ftype == "mRNA":
                tmpcounter += 1
                tmpnewgeneid = "%s_gene_%i" % (annoversion, tmpcounter)
                tmpnewgene = Feature(seqid=feature.seqid, source=feature.source, ftype="gene", start=feature.start, end=feature.end, score=feature.score, strand=feature.strand, phase=feature.phase)
                tmpnewgene.identifier = tmpnewgeneid
                tmpnewgene.name = tmpnewgeneid
                tmpnewgene.alias = feature.parent.alias
                tmpnewgene.notes = feature.parent.notes
                tmpnewgene.target = feature.parent.target
                tmpnewgene.primary_confidence_class = feature.parent.primary_confidence_class
                tmpnewgene.secondary_condidence_class = feature.parent.secondary_condidence_class
                feature.parent = tmpnewgene
                tmpnewgene.features = [feature]
                tmpnewgenes.append(tmpnewgene)
        
        #2) merge genes with overlapping CDS (features need to be sorted)
        opencdss = {}
        removegeneids = set([])
        for feature in self.features:
            if feature.ftype=="CDS":
                tmpopencdss = []
                opengeneid = "none"
                currentgene = feature.parent.parent
                # if there are no open CDS for current seqid initialze empty array
                if not feature.seqid in opencdss.keys():
                    opencdss[feature.seqid] = []
                # go through all open cds and keep if still open; set new gene to last open cds (gene are same for all open CDS on same strand)
                for opencds in opencdss[feature.seqid]:
                    if opencds.end>=feature.start:
                        tmpopencdss.append(opencds)
                        if opencds.strand==feature.strand:
                            opengene = opencds.parent.parent
                            opengeneid = opengene.identifier
                # set new gene to last open cds gene
                if currentgene.identifier!=opengeneid and opengeneid!="none":
                    tmpstart=math.inf
                    tmpend=0
                    for tmpmrna in currentgene.features:
                        tmpmrna.parent = opengene
                        opengene.features.append(tmpmrna)
                        tmpstart = min(tmpstart, tmpmrna.start)
                        tmpend = max(tmpend, tmpmrna.end)
                    opengene.start = min(tmpstart, opengene.start)
                    opengene.end = max(tmpend, opengene.end)
                    if currentgene.source!=opengene.source:
                        opengene.source = "multiple"
                    currentgene.mrnas = []
                    removegeneids.add(currentgene.identifier)
                tmpopencdss.append(feature)
                opencdss[feature.seqid] = tmpopencdss
        
        #3) update features and return object
        newfeatures = []
        newgenes = {}
        newseqids = {}
        for feature in self.features:
            if feature.ftype!="gene":
                newfeatures.append(feature)
        for gene in tmpnewgenes:
            if not gene.identifier in removegeneids:
                if not gene.seqid in newseqids.keys():
                    newseqids[gene.seqid] = []
                newgenes[gene.identifier] = gene
                newseqids[gene.seqid].append(gene)
                newfeatures.append(gene)
        self.genes = newgenes
        self.seqids = newseqids
        self.features = sorted(newfeatures)
        return(self)
    
    def collapseMrnas(self):
        """
        This function removes redundant mRNAs
        """
        
        newfeatures = []
        newmrnas = {}
        # go through all genes
        for geneid in self.genes:
            gene = self.genes[geneid]
            newfeatures.append(gene)
            tmp_keeptranscripts = []
            # go through all mRNAs
            for mrna1 in gene.features:
                isequal = False
                # check if there is already equal mRNA in set of new mRNAs
                for mrna2 in tmp_keeptranscripts:
                    if mrna1 == mrna2:
                        isequal = True
                if not isequal:
                    tmp_keeptranscripts.append(mrna1)
            # set new mRNAs
            gene.features = tmp_keeptranscripts
            # add features to newfeatures (those to keep)
            for mrna in tmp_keeptranscripts:
                newfeatures.append(mrna)
                newmrnas[mrna.identifier] = mrna
                newfeatures += mrna.features
        self.features = sorted(newfeatures)
        self.mrnas = newmrnas
        return(self)
    
    def writeGff3Genes(self, path):
        with open(path, "w") as outgff:
            writer = csv.writer(outgff, delimiter="\t", quotechar="#", quoting = csv.QUOTE_NONE)
            for feature in self.features:
                writer.writerow(feature.getLine())
    
    def printGeneStats(self, path):
        with open(path, "w") as outfile:
            rowpattern = {"id":"none", "source":"none", "seqid":"none", "start":0, "end":0, "ntranscripts":0, "primconf":""}
            variables = ["id", "source", "seqid", "start", "end", "ntranscripts", "primconf"]
            writer = csv.DictWriter(outfile, fieldnames=variables)
            writer.writeheader()
            for geneid in self.genes:
                gene = self.genes[geneid]
                outdata = dict(rowpattern)
                outdata["id"] = geneid
                outdata["source"] = gene.source
                outdata["seqid"] = gene.seqid
                outdata["start"] = gene.start
                outdata["end"] = gene.end
                outdata["ntranscripts"] = len(gene.features)
                outdata["primconf"] = gene.primary_confidence_class
                writer.writerow(outdata)
    
    def printTranscriptsStats(self, path, includetargets=False):
        with open(path, "w") as outfile:
            rowpattern = {"id":"none", "gene": "none", "source":"none", "seqid":"none", "start":0, "end":0, "bpcdss":0, "ncdss":0, "primconf":"", "secconf":""}
            variables = ["id", "gene", "source", "seqid", "start", "end", "bpcdss", "ncdss", "primconf", "secconf"]
            if includetargets:
                rowpattern["target"] = ""
                variables.append("target")
            writer = csv.DictWriter(outfile, fieldnames=variables)
            writer.writeheader()
            for mrnaid in self.mrnas:
                mrna = self.mrnas[mrnaid]
                outdata = dict(rowpattern)
                outdata["id"] = mrnaid
                outdata["gene"] = mrna.parent.identifier
                outdata["source"] = mrna.source
                outdata["seqid"] = mrna.seqid
                outdata["start"] = mrna.start
                outdata["end"] = mrna.end
                outdata["primconf"] = mrna.primary_confidence_class
                outdata["secconf"] = mrna.secondary_condidence_class
                tmpbpcdss = 0
                tmpncdss = 0
                for cds in mrna.features:
                    if cds.ftype=="CDS":
                        tmpncdss += 1
                        tmpbpcdss += (cds.end-cds.start)+1
                outdata["ncdss"] = tmpncdss
                outdata["bpcdss"] = tmpbpcdss
                if includetargets:
                    outdata["target"] = mrna.target
                writer.writerow(outdata)
    
    def getHcGff3Genes(self):
        newfeatures = []
        newgenes = {}
        newseqids = {}
        newmrnas = {}
        for geneid in self.genes:
            gene = self.genes[geneid]
            if gene.primary_confidence_class=="HC":
                newfeatures.append(gene)
                newgenes[gene.identifier] = gene
                if not gene.seqid in newseqids:
                    newseqids[gene.seqid] = []
                newseqids[gene.seqid].append(gene)
                for mrna in gene.features:
                    newmrnas[mrna.identifier] = mrna
                    newfeatures.append(mrna)
                    newfeatures += mrna.features
        newanno = GeneAnnotation()
        newanno.features = sorted(newfeatures)
        newanno.genes = newgenes
        newanno.seqids = newseqids
        newanno.mrnas = newmrnas
        return newanno
    
    def getLcGff3Genes(self):
        newfeatures = []
        newgenes = {}
        newseqids = {}
        newmrnas = {}
        for geneid in self.genes:
            gene = self.genes[geneid]
            if gene.primary_confidence_class=="LC":
                newfeatures.append(gene)
                newgenes[gene.identifier] = gene
                if not gene.seqid in newseqids:
                    newseqids[gene.seqid] = []
                newseqids[gene.seqid].append(gene)
                for mrna in gene.features:
                    newmrnas[mrna.identifier] = mrna
                    newfeatures.append(mrna)
                    newfeatures += mrna.features
        newanno = GeneAnnotation()
        newanno.features = sorted(newfeatures)
        newanno.genes = newgenes
        newanno.seqids = newseqids
        newanno.mrnas = newmrnas
        return newanno

