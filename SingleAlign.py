#!/usr/bin/env python

import re
import pprint

import __Utils
import SingleSamAlign

##### Auxilary functions:
#

compSeq = __Utils.compSeq
bioExplainCigar = __Utils.bioExplainCigar

#
#####

class SingleAlign:
	__slots__ = ["primary", "supple"]
	
	def __init__(self, ref):
		self.primary = ref["primary"].toSingleAlignRef()
		self.supple  = ref["supple"].toSingleAlignRef()


	def primaryAlignment(self):
		return self.primary


	def suppleAlignment(self):
		return self.supple


	def para(self):
		return self.primary["para"]


	def getAlign(self, id):
		if id == "primary":
			return self.primary

		if id == "supple":
			return self.supple

		return None


	def rnam(self, id):
		if id == "primary":
			return self.primary["rnam"]

		if id == "supple":
			return self.supple["rnam"]

		return None


	def flag(self, id):
		if id == "primary":
			return self.primary["flag"]

		if id == "supple":
			return self.supple["flag"]

		return None


	def chro(self, id):
		if id == "primary":
			return self.primary["chro"]

		if id == "supple":
			return self.supple["chro"]

		return None


	def strand(self, id, cha = None):
		if id == "primary":
			if cha != None:
				self.primary["strand"] = cha

			return self.primary["strand"]

		if id == "supple":
			if cha != None:
				self.supple["strand"] = cha

			return self.supple["strand"]

		return None


	def human_strand(self):
		chrVirus = self.para().virus

		human = None
		if self.chro("primary") == chrVirus:
			human = "supple"

		if self.chro("supple") == chrVirus:
			human = "primary"

		return self.strand(human)


	def start(self, id, sta = None):
		if id == "primary":
			if sta != None:
				self.primary["start"] = sta

			return self.primary["start"]

		if id == "supple":
			if sta != None:
				self.supple["start"] = sta

			return self.supple["start"]

		return None


	def stop(self, id, sto = None):
		if id == "primary":
			if sto != None:
				self.primary["stop"] = sto

			return self.primary["stop"]

		if id == "supple":
			if sto != None:
				self.supple["stop"] = sto

			return self.supple["stop"]

		return None


	def mapq(self, id):
		if id == "primary":
			return self.primary["mapq"]

		if id == "supple":
			return self.supple["mapq"]

		return None


	def cigar(self, id, cigar = None):
		if id == "primary":
			if cigar != None:
				self.primary["cigar"] = cigar

			return self.primary["cigar"]

		if id == "supple":
			if cigar != None:
				self.supple["cigar"] = cigar

			return self.supple["cigar"]

		return None


	def formatCigar(self, id):
		if id != "primary" and id != "supple":
			return None

		ss = "self.cigar(\"%s\")" %(id)
		cigar = eval(ss)

		pat = re.compile("[0-9]+[SH][0-9]+M[0-9]+[SH]")
		if re.search(pat, ss) != None:
			return cigar

		pat = re.compile("[0-9]+[SH][0-9]+M")
		if re.search(pat, cigar) != None:
			return "%s0S" %(cigar)

		pat = re.compile("[0-9]+M[0-9]+[SH]")
		if re.search(pat, cigar) != None:
			return "0S%s" %(cigar)


	def sequ(self, id, ss = None):
		if id == "primary":
			if ss != None:
				self.primary["sequ"] = ss

			return self.primary["sequ"]

		if id == "supple":
			if ss != None:
				self.supple["sequ"] = ss

			return self.supple["sequ"]

		return None


	def qual(self, id, ss = None):
		if id == "primary":
			if ss != None:
				self.primary["qual"] = ss

			return self.primary["qual"]

		if id == "supple":
			if ss != None:
				self.supple["qual"] = ss

			return self.supple["qual"]

		return None


	# def revQual(self, id):
	# 	# ss = f"self.qual(\"%s\")"
	# 	ss = "self.qual('%s')" %(id)
	# 	qual = eval(ss)
	# 	qual = qual[::-1]

	# 	# ss = f"self.qual({id}, {qual})"
	# 	ss = "self.qual('%s', '%s')" %(id, qual)
	# 	eval(ss)

	# 	return 0

	def revQual(self, id):
		qual = self.qual(id)
		self.qual(id, qual[::-1])

		return 0


	# def revSequQual(self, id):
	# 	ss = "self.sequ(\"%s\")" %(id)
	# 	sequ = eval(ss)[::-1]

	# 	ss = "self.sequ(\"%s\", \"%s\")" %(id, sequ)
	# 	eval(ss)

	# 	ss = "self.revQual(\"%s\")" %(id)
	# 	eval(ss)

	# 	return 0

	def revSequQual(self, id):
		sequ = self.sequ(id)
		self.sequ(id, sequ[::-1])

		self.revQual(id)

		return 0


	# def compSequ(self, id):
	# 	ss = "self.sequ(\"%s\")" %(id)
	# 	seq = eval(ss)

	# 	ss = compSeq(seq)
	# 	eval("self.sequ(\"%s\", \"%s\")" %(id, ss))


	def compSequ(self, id):
		seq = self.sequ(id)
		ss  = compSeq(seq)
		self.sequ(id, ss)


	def set_reverse_strand(self, id):
		cha = self.strand(id)

		if(cha == "+"):
			self.strand(id, "-")

		if(cha == "-"):
			self.strand(id, "+")

		self.compSequ(id)
		return 0


	def has_supplement(self):
		if self.rnam("supple") != "NA":
			return 1
		else:
			return 0


	def is_hs_virus(self):
		chrVirus = self.para().virus

		chrp = self.chro("primary")
		chrs = self.chro("supple")

		if chrp != chrVirus and chrs != chrVirus:
			return 1

		if chrp == chrVirus and chrs == chrVirus:
			return 1

		return 0


	def bioExplainCigar(self, id):
		ss = "self.cigar(\"%s\")" %(id)
		cigar_str = eval(ss)
		return bioExplainCigar(cigar_str)


	def reverseCigar(self, id):
		cigar_ref = self.bioExplainCigar(id)

		keys = list(cigar_ref.keys())
		keys_rev = keys[::-1]

		res = {}
		for k, kv in zip(keys, keys_rev):
			res[k] = cigar_ref[kv]

		return res


	def set_reverse_cord(self, id):
		sta = self.start(id)
		sto = self.stop(id)

		self.start(id, sto)
		self.stop(id, sta)


	def set_reverse_cigar(self, id):
		cigar_rev = self.reverseCigar(id)

		res = ""
		for key_tmp in cigar_rev.keys():
			nbase  = cigar_rev[key_tmp]["nbase"]
			symbol = cigar_rev[key_tmp]["symbol"]

			res += str(nbase) + symbol

		self.cigar(id, res)
		self.set_reverse_cord(id)
		self.revSequQual(id)


	def figure_cigar(self):
		chrVirus = self.para().virus

		cigar_p = self.formatCigar("primary")
		cigar_s = self.formatCigar("supple")

		np = None
		ns = None

		res = {
			"left"       : "",
			"left_id"    : "",
			"left_rev"   : "",
			"left_spec"  : "",
			"right"      : "",
			"right_id"   : "",
			"right_rev"  : "",
			"right_spec" : "",
		}

		pat = re.compile("[0-9]+[SH]")

		resu = re.match(pat, cigar_p)
		if resu == None:
			print("Can't figure out np!")
			np = 0
		else:
			np = int(resu.group()[:-1])

		resu = re.match(pat, cigar_s)
		if resu == None:
			print("Can't figure out ns!")
			ns = 0
		else:
			ns = int(resu.group()[:-1])

		if np == 0 and ns == 0:
			print("Both np and ns are 0")
			return "NA"

		if np >= ns:
			res["left"] = self.suppleAlignment()
			res["left_id"] = "supple"
			res["right"] = self.primaryAlignment()
			res["right_id"] = "primary"

			if self.start("primary") < self.stop("primary"):
				res["right_rev"] = 0
			else:
				res["right_rev"] = 1

			if self.start("supple") < self.stop("supple"):
				res["left_rev"] = 0
			else:
				res["left_rev"] = 1

			if self.chro("primary") == chrVirus:
				res["left_spec"] = "human"
				res["right_spec"] = "virus"
			else:
				res["left_spec"] = "virus"
				res["right_spec"] = "human"
		else:
			res["left"] = self.primaryAlignment()
			res["left_id"] = "primary"
			res["right"] = self.suppleAlignment()
			res["right_id"] = "supple"

			if self.start("primary") < self.stop("primary"):
				res["left_rev"] = 0
			else:
				res["left_rev"] = 1

			if self.start("supple") < self.stop("supple"):
				res["right_rev"] = 0
			else:
				res["right_rev"] = 1

			if self.chro("primary") == chrVirus:
				res["left_spec"] = "virus"
				res["right_spec"] = "human"
			else:
				res["left_spec"] = "human"
				res["right_spec"] = "virus"

		cigar_ref = {
			"primary" : cigar_p,
			"supple"  : cigar_s,
		}

		return res, cigar_ref


	def alignConfiguration(self):
		chrVirus = self.para().virus

		if self.chro("supple") == "NA":
			return "NA"

		align_human = None
		align_virus = None
		if self.chro("primary") == chrVirus:
			align_human = "supple"
			align_virus = "primary"
		else:
			align_human = "primary"
			align_virus = "supple"

		if self.strand(align_virus) != self.strand(align_human):
			self.set_reverse_cigar(align_virus)

		config, cigar_ref = self.figure_cigar()
		cigar_left  = cigar_ref[ config["left_id"]  ]
		cigar_right = cigar_ref[ config["right_id"] ]

		homolog = {
			"length" : 0,
			"human"  : {
				"chro"   : "NA",
				"strand" : ".",
				"start"  : 0,
				"stop"   : 0,
			},
			"virus"  : {
				"chro"   : "NA",
				"strand" : ".",
				"start"  : 0,
				"stop"   : 0,
			},
		}

		sl_left  = None
		ma_left  = None
		sl_right = None

		pat = re.compile("[0-9]+[SH]")
		sl_left = int(re.match(pat, cigar_left).group()[:-1])

		pat = re.compile("[0-9]+M")
		ma_left = int(re.search(pat, cigar_left).group()[:-1])

		pat = re.compile("[0-9]+[SH]")
		sl_right = int(re.match(pat, cigar_right).group()[:-1])

		hlen = sl_left + ma_left - sl_right
		homolog["length"] = hlen

		if hlen == 0:
			config["homolog"] = homolog
			return config

		strand_flag = 0
		if config["left_rev"] == 1 or config["right_rev"] == 1:
			strand_flag = 1

		if strand_flag == 0:
			sta_left  = config["left"]["start"]
			sta_right = config["right"]["start"]

			staL, stoL = sta_left + ma_left - hlen, sta_left + ma_left - 1
			staR, stoR = sta_right, sta_right + hlen - 1

			human = None
			virus = None
			if config["left_spec"] == "human":
				homolog["human"]["start"] = staL
				homolog["human"]["stop"]  = stoL
				homolog["virus"]["start"] = staR
				homolog["virus"]["stop"]  = stoR

				human = "left"
				virus = "right"
			else:
				homolog["human"]["start"] = staR
				homolog["human"]["stop"]  = stoR
				homolog["virus"]["start"] = staL
				homolog["virus"]["stop"]  = stoL

				human = "right"
				virus = "left"

			homolog["human"]["chro"] = config[human]["chro"]
			homolog["human"]["strand"] = config[human]["strand"]
			homolog["virus"]["chro"] = config[virus]["chro"]
			homolog["virus"]["strand"] = config[virus]["strand"]
		else:
			if config["left_spec"] == "human":
				sto_human = config["left"]["stop"]
				sta_virus = config["right"]["start"]

				staHuman, stoHuman = sto_human - hlen + 1, sto_human
				staVirus, stoVirus = sta_virus, sta_virus - hlen + 1

				homolog["human"]["chro"] = config["left"]["chro"]
				homolog["human"]["strand"] = config["left"]["strand"]
				homolog["human"]["start"] = staHuman
				homolog["human"]["stop"] = stoHuman

				homolog["virus"]["chro"] = config["right"]["chro"]
				homolog["virus"]["strand"] = config["right"]["strand"]
				homolog["virus"]["start"] = staVirus
				homolog["virus"]["stop"] = stoVirus
			else:
				sta_human = config["right"]["start"]
				sto_virus = config["left"]["stop"]

				staHuman, stoHuman = sta_human, sta_human + hlen - 1
				staVirus, stoVirus = sto_virus + hlen - 1, sto_virus

				homolog["human"]["chro"] = config["right"]["chro"]
				homolog["human"]["strand"] = config["right"]["strand"]
				homolog["human"]["start"] = staHuman
				homolog["human"]["stop"] = stoHuman

				homolog["virus"]["chro"] = config["left"]["chro"]
				homolog["virus"]["strand"] = config["left"]["strand"]
				homolog["virus"]["start"] = staVirus
				homolog["virus"]["stop"] = stoVirus

		config["homolog"] = homolog
		return config


if __name__ == "__main__":
	import __ParseArgs
	import SingleSamAlign

	args = __ParseArgs.parseArgs()

	title = {
		"C1" : "Case_1: Human (+) / HBV (+);",
		"C2" : "Case_2: Human (+) / HBV (-);",
		"C3" : "Case_3: Human (-) / HBV (+);",
		"C4" : "Case_4: Human (-) / HBV (-);",
		"C5" : "Case_5: Human (+) / HBV (+);",
		"C6" : "Case_6: Human (+) / HBV (-);",
		"C7" : "Case_7: Human (-) / HBV (+);",
		"C8" : "Case_8: Human (-) / HBV (-);",
	}

	rec = {}
	reads = []

	with open("Test_Data_SingleAlign.sam", "r") as f:
		for line in f:
			pat = re.compile("_C[1-8]#\t")

			if re.search(pat, line) != None:
				continue

			tmp = line.split("\t")

			if not tmp[0] in rec:
				reads.append(tmp[0])
				rec[ tmp[0] ] = []

			rec[ tmp[0] ].append(line[:-1])

	for read_tmp in reads:
		id = read_tmp[-2:]
		print(">%s" %(title[id]))

		sam1 = SingleSamAlign.SingleSamAlign(rec[read_tmp][0], args)
		sam2 = SingleSamAlign.SingleSamAlign(rec[read_tmp][1], args)

		print("\n ----- Read 1 -----")
		sam1.printClass()

		print("\n ----- Read 2 -----")
		sam2.printClass()

		ref = {
			sam1.getAlignID() : sam1,
			sam2.getAlignID() : sam2,
		}
		align = SingleAlign(ref)

		print("\n ----- primary alignment -----")
		alignp = align.primaryAlignment()
		print(alignp)

		print("\n ----- supplementary alignment -----")
		aligns = align.suppleAlignment()
		print(aligns)

		print("\n ----- parameters -----")
		args = align.para()
		print(args)

		print("\n ----- test for getAlign -----")
		testp = align.getAlign("primary")
		tests = align.getAlign("supple")
		print(" ----- test OK -----")

		print("The read name is: %s" %(align.rnam("primary")))
		print("Human strand is: %s" %(align.human_strand()))

		print("The primary CIGAR is: %s" %(align.cigar("primary")))
		print("The formated CIGAR is: %s" %(align.formatCigar("primary")))

		print("The supplementary CIGAR is: %s" %(align.cigar("supple")))
		print("The formated CIGAR is: %s" %(align.formatCigar("supple")))

		# print("The primary sequence is: %s" %(align.sequ("primary")))
		# print("The primary quality is: %s" %(align.qual("primary")))

		# align.revSequQual("primary")
		# print("The primary sequence is: %s" %(align.sequ("primary")))
		# print("The primary quality is: %s" %(align.qual("primary")))

		# print("The primary sequence is: %s" %(align.sequ("primary")))
		# align.compSequ("primary")
		# print("The primary sequence is: %s" %(align.sequ("primary")))

		# print("The strand for primary is: %s" %(align.strand("primary")))
		# align.set_reverse_strand("primary")
		# print("The reverse strand for primary is: %s" %(align.strand("primary")))

		# if align.has_supplement():
		# 	print("There's supplementary mapping.")
		
		# if not align.is_hs_virus():
		# 	print("The mapping is chimeric.")
		
		# cigar_ref = align.bioExplainCigar("primary")
		# print(cigar_ref)

		# print("The primary CIGAR is: %s" %(align.cigar("primary")))
		# cigar_rev = align.reverseCigar("primary")
		# print(cigar_rev)

		# print("The start of primary mapping: %d" %(align.start("primary")))
		# print("The stop of primary mapping: %d" %(align.stop("primary")))
		# align.set_reverse_cord("primary")
		# print("The start of primary mapping: %d" %(align.start("primary")))
		# print("The stop of primary mapping: %d" %(align.stop("primary")))

		# align.set_reverse_cigar("primary")
		# print("The reverse cigar is: %s" %(align.cigar("primary")))

		# config, cigar_ref = align.figure_cigar()
		# print(config)
		# print(cigar_ref)

		config = align.alignConfiguration()
		pp = pprint.PrettyPrinter(indent = 4)
		pp.pprint(config)

		print("\n -------------------- \n")

