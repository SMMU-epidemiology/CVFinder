#!/usr/bin/env python

import re
import cigar
# import __Utils

class SingleSamAlign:
	__slots__ = [
		"rnam", "flag", "chro",
		"cord", "mapq", "cigar",
		"sequ", "qual", "args",
	]
	
	def __init__(self, rec, args):
		if rec == "":
			self.rnam  = "NA"
			self.flag  = 0
			self.chro  = "NA"
			self.cord  = 0
			self.mapq  = 60
			self.cigar = "NA"
			self.sequ  = "NA"
			self.qual  = "NA"
			self.args  = -1
		else:
			tmp = rec.split("\t")

			self.rnam = tmp[0]
			self.flag = int(tmp[1])

			chro = tmp[2]
			if not chro.startswith("chr"):
				chro = "chr" + chro

			self.chro  = chro
			self.cord  = int(tmp[3])
			self.mapq  = int(tmp[4])
			self.cigar = tmp[5]
			self.sequ  = tmp[9].upper()
			self.qual  = tmp[10]
			self.args  = args
	

	def printClass(self):
		print("Read name is: %s" %(self.rnam))
		print("The flag is: %d" %(self.flag))
		print("The chromosome is: %s" %(self.chro))
		print("The coordinate is: %d" %(self.cord))
		print("The mapping quality is: %d" %(self.mapq))
		print("The CIGAR is: %s" %(self.cigar))
		print("The sequence is: %s" %(self.sequ))
		print("The sequencing quality is: %s" %(self.qual))
		print("The arguments are:")
		print(self.args)


	def strand(self):
		if self.flag & 16:
			return "-"
		else:
			return "+"


	def revQual(self):
		qual = self.qual
		self.qual = qual[::-1]
	

	def revSequQual(self):
		sequ = self.sequ
		self.sequ = sequ[::-1]
		self.revQual()


	def getSpec(self):
		if self.chro == "NA":
			return "NA"

		if self.chro == self.args.virus:
			return "Virus"
		else:
			return "Human"


	def getAlignID(self):
		if self.flag & 2048:
			return "supple"
		else:
			return "primary"


	def is_dummy(self):
		if self.rnam == "NA":
			return 1
		else:
			return 0


	def is_not_simple(self):
		pat = re.compile("D|I")
		if re.search(pat, self.cigar) == None:
			return 0
		else:
			return 1


	def idName(self):
		id = "NA"
		sa = "NA"

		if not self.flag & 1:
			id = 1

		if self.flag & 64:
			id = 1

		if self.flag & 128:
			id = 2

		if self.flag & 2048:
			sa = "_sa"
		else:
			sa = ""

		return "read%s%s" %(id, sa)


	def is_OK(self):
		"""
		Test if the flag has: 4 (read unmapped), 256 (not primary alignment),
		512 (read fails platform/vendor quality checks) and 1024 (read is PCR or optical duplicate);
		"""
		# if self.flag & 1796:
		if self.flag & 1792:
			return 0
		else:
			return 1


	def isSingleEndSeq(self):
		if self.flag & 1:
			return 0
		else:
			return 1


	def bioExplainCigar(self):
		ref = cigar.Cigar(self.cigar)
		res = {}

		for i, ref_tmp in enumerate(ref.items()):
			res[i] = {
				"nbase" : ref_tmp[0],
				"symbol": ref_tmp[1],
			}

		return res
	

	def fixCigar(self):
		"""
		Fix Alignments with Deletion / Insertion
		1. Remove Inertions from sequence and cigar;
		2. Add dashes (i.e. "-") for deletions on the sequence;
		3. Objects are checked by "is_not_simple" first.
		4. Check: OK
		"""
		cigar_ref = self.bioExplainCigar()

		nseg = 0
		nins = "NA"
		ndel = "NA"

		for key_tmp in cigar_ref.keys():
			syn = cigar_ref[key_tmp]["symbol"]

			if syn == "M" or syn == "S":
				nseg += cigar_ref[key_tmp]["nbase"]
				continue

			if syn == "I":
				nins = cigar_ref[key_tmp]["nbase"]
				seg1 = self.sequ[0:nseg]
				seg2 = self.sequ[nseg + nins:]
				self.sequ = seg1 + seg2

				quaf = self.qual[0:nseg]
				quas = self.qual[nseg + nins:]
				self.qual = quaf + quas

				continue

			if syn == "D":
				ndel = cigar_ref[key_tmp]["nbase"]

				seg1 = self.sequ[0:nseg]
				seg2 = self.sequ[nseg:]
				self.sequ = seg1 + "-"*ndel + seg2

				quaf = self.qual[0:nseg]
				quas = self.qual[nseg:]
				self.qual = quaf + "-"*ndel + quas

				# After adding dashes, deletion is equivalent to match. 
				nseg += ndel

				continue

		# Fix cigar string
		# Head
		#
		cigar_head = "NA"
		nhead = 0
		chead = ""
		key_head = list(cigar_ref.keys())[0]

		if cigar_ref[key_head]["symbol"] == "H" or cigar_ref[key_head]["symbol"] == "S":
			nhead = cigar_ref[key_head]["nbase"]
			chead = cigar_ref[key_head]["symbol"]

		# Tail
		#
		cigar_tail = "NA"
		ntail = 0
		ctail = ""
		key_tail = list(cigar_ref.keys())[-1]

		if cigar_ref[key_tail]["symbol"] == "H" or cigar_ref[key_tail]["symbol"] == "S":
			ntail = cigar_ref[key_tail]["nbase"]
			ctail = cigar_ref[key_tail]["symbol"]

		# Middle
		#
		nseq = len(self.sequ)
		if chead == "S":
			nseq -= nhead

		if ctail == "S":
			nseq -= ntail

		cigar_fix = ""
		if nhead > 0:
			cigar_fix += str(nhead) + chead

		cigar_fix += str(nseq) + "M"

		if ntail > 0:
			cigar_fix += str(ntail) + ctail

		self.cigar = cigar_fix


	def getMatchLength(self):
		cigar_ref = self.bioExplainCigar()
		res = 0

		for key_tmp in cigar_ref.keys():
			if cigar_ref[key_tmp]["symbol"] == "M":
				res += cigar_ref[key_tmp]["nbase"]

		return res


	def readLength(self):
		cigar_str = cigar.Cigar(self.cigar)
		return len(cigar_str)


	def tailToHead(self, obje):
		res = {
			"stat" : 0,
			"merge": "NA",
		}

		vec = {
			"tail": "NA",
			"head": "NA",
		}

		############### Some Simple Test:

		chra = self.chro
		chrb = obje.chro

		chrVirus = self.args.virus
		if chra != chrVirus or chrb != chrVirus:
			return res

		strand1 = self.strand()
		strand2 = obje.strand()

		if strand1 != strand2:
			return res

		###############

		sta1 = self.cord
		mat1 = self.getMatchLength()
		sta2 = obje.cord
		mat2 = obje.getMatchLength()

		if (not sta1 != 1) and (not sta2 != 1):
			return res

		nbas = max(self.readLength(), obje.readLength())

		cigar1 = self.cigar
		cigar2 = obje.cigar

		if cigar1.find("H") == -1:
			nbas = len(self.sequ)
		elif cigar2.find("H") == -1:
			nbas = len(obje.sequ)
		else:
			nbas = len(self.sequ) + len(obje.sequ)

		###########   A Rigid Standard:
		if mat1 + mat2 != nbas:
			return res

		###########

		if sta1 == 1:
			if sta2 + mat2 -1 != self.args.vlen:
				return res

			res["stat"] = 1
			vec["tail"] = obje
			vec["head"] = self
		elif sta2 == 1:
			if sta1 + mat1 - 1 != self.args.vlen:
				return res

			res["stat"] = 1
			vec["tail"] = self
			vec["head"] = obje
		else:
			return res

		# ref = {
		# 	"rnam" : self.rnam,
		# 	"flag" : 0,
		# 	"chro" : chrVirus
		# 	"cord" : 0,
		# 	"mapq" : 60,
		# 	"cigar": str(nbas) + "M",
		# 	"sequ" : "NA",
		# 	"qual" : "NA",
		# 	"vlen" : self.args.virusLength,
		# }

		ref = SingleSamAlign("", self.args)

		ref.rnam  = self.rnam
		ref.chro  = self.chro
		ref.flag = self.flag & obje.flag
		sta = vec["tail"].cord
		ref.cord = -sta
		ref.mapq = max(self.mapq, obje.mapq)

		cigar1 = self.cigar
		cigar2 = obje.cigar

		if cigar1.find("H") == -1:
			ref.sequ = self.sequ
			ref.qual = self.qual
		elif cigar2.find("H") == -1:
			ref.sequ = obje.sequ
			ref.qual = obje.qual
		else:
			ref.sequ = vec["tail"].sequ + vec["head"].sequ
			ref.qual = vec["tail"].qual + vec["head"].qual

		nbas = len(ref.sequ)
		ref.cigar = str(nbas) + "M"

		ref.args = self.args
		res["merge"] = ref

		return res


	def headOrTail(self):
		if self.chro != self.args.virus:
			return -1

		if abs(self.cord)-1 <= self.args.vlen-abs(self.cord):
			return 0
		else:
			return 1


	def toSingleAlignRef(self):
		res = {
			"rnam"  : self.rnam,
			"flag"  : self.flag,
			"chro"  : self.chro,
			"strand": ".",
			"start" : self.cord,
			"stop"  : 0,
			"mapq"  : self.mapq,
			"cigar" : self.cigar,
			"sequ"  : self.sequ,
			"qual"  : self.qual,
			"circ"  : 0,
			"para"  : self.args,
		}

		if self.is_dummy():
			return res

		# if self.flag & 16:
		# 	res["strand"] = "-"
		# else:
		# 	res["strand"] = "+"

		res["strand"] = self.strand()

		if self.cord > 0:
			res["stop"] = res["start"] + self.getMatchLength() - 1

		if self.cord < 0:
			ntail = self.args.vlen + self.cord + 1
			nhead = self.readLength() - ntail
			res["stop"] = nhead
			res["circ"] = 1

		return res


if __name__ == "__main__":
	import __ParseArgs
	args = __ParseArgs.parseArgs()

	# with open("Test_Data_SingleSamAlign.sam", "r") as f:
	with open("Test_Data.sam", "r") as f:
		while(True):
			rec = f.readline().strip()

			if rec.startswith("@"):
				continue
			else:
				break

		# rec = f.readline().strip()
		test = SingleSamAlign(rec, args)
		test.printClass()

		# test.revSequQual()
		print("Test for revSequQual:")
		print("Now the sequence is: %s" %(test.sequ))
		print("Now the sequencing quality is: %s" %(test.qual))
		print("The Specie is: %s" %(test.getSpec()))
		print("The read is: %s" %(test.getAlignID()))
		print("The read name is: %s" %(test.idName()))
		print("The read is OK: %d" %(test.is_OK()))
		print("The read is single-ended: %d" %(test.isSingleEndSeq()))
		print("The cigar dict is:")
		print(test.bioExplainCigar())

		print("\n -------------------------------\n")
		print("Before fixing the CIGAR flag:")
		# rec = f.readline().strip()
		# temp = SingleSamAlign(rec, args)
		# temp.printClass()
		test.printClass()

		print("\n -------------------------------\n")
		print("After fixing the CIGAR flag:")
		test.fixCigar()
		test.printClass()

		# print("\n -------------------------------\n")
		# print("Before fixing the CIGAR flag:")
		# rec = f.readline().strip()
		# temp = SingleSamAlign(rec, args)
		# temp.printClass()
		# print("The match length is: %d" %(temp.getMatchLength()))
		# print("The read length is: %d" %(temp.readLength()))

		# print("\n -------------------------------\n")
		# print("After fixing the CIGAR flag:")
		# temp.fixCigar()
		# temp.printClass()

		exit(0)

		print("\n -------------------------------\n")
		print("Test for <tailToHead>:")
		rec = f.readline().strip()
		read1 = SingleSamAlign(rec, args)
		rec = f.readline().strip()
		read2 = SingleSamAlign(rec, args)

		read1.printClass()
		print("Read 1 <headOrTail>: ", read1.headOrTail())

		read2.printClass()
		print("Read 2 <headOrTail>: ", read2.headOrTail())

		print("\n -------------------------------\n")
		print("After mergeing:")
		res = read1.tailToHead(read2)
		print(res)
		res["merge"].printClass()

		print("\n -------------------------------\n")
		print("Test for <toSingleAlignRef>:")
		ref = res["merge"].toSingleAlignRef()
		print(ref)
