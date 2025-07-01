#!/usr/bin/env python

import FragmentAlign
import FragmentCluster

import __Utils


##### Auxilary functions:
#

getChromosomes = __Utils.getChromosomes

#
#####


class IdeogramAlign:
	__slots__ = [
		"orientation", "data",
	]
	
	def __init__(self, para):
		chrs_ref = getChromosomes(para)
		chrs     = chrs_ref["human"]

		config = [
			"Human:+ | Virus:+", "Human:+ | Virus:-",
			"Virus:+ | Human:+", "Virus:- | Human:+"
		]
		self.orientation = config

		ref = {
			"Human:+ | Virus:+" : {},
			"Human:+ | Virus:-" : {},
			"Virus:+ | Human:+" : {},
			"Virus:- | Human:+" : {},
		}


		for chro in chrs:
			ref["Human:+ | Virus:+"][chro] = dict(
				frags  = [],
				nfrags = 0,
			)

			ref["Human:+ | Virus:-"][chro] = dict(
				frags  = [],
				nfrags = 0,
			)

			ref["Virus:+ | Human:+"][chro] = dict(
				frags  = [],
				nfrags = 0,
			)

			ref["Virus:- | Human:+"][chro] = dict(
				frags  = [],
				nfrags = 0,
			)

		self.data = ref

	
	def getAllOrientation(self):
		return self.orientation


	def getAllChromosome(self):
		res = self.data["Human:+ | Virus:+"].keys()
		return list(res)


	def frags(self, orientation, chro):
		res = self.data[orientation][chro]["frags"]
		return res


	def nfrags(self, orientation, chro):
		res = self.data[orientation][chro]["nfrags"]
		return res


	def addFragment(self, frag):
		if not isinstance(frag, FragmentAlign.FragmentAlign):
			return

		orientation = frag.orientation
		config      = frag.configPosi()
		chro        = frag.chro(config["human"])
		idx         = self.nfrags(orientation, chro)

		self.data[orientation][chro]["frags"].append(frag)
		self.data[orientation][chro]["nfrags"] = idx + 1


	def sortFrags(self, orient, chro, spec = None):
		if spec == None:
			spec = "human"

		frags = self.frags(orient, chro)
		if self.nfrags(orient, chro) < 2:
			return frags

		# posi = frags[0].configPosi()[spec]
		frags_sort = sorted(
			frags,
			key = lambda x: x.figureBreakpoint(spec)["cord"]
		)
		return frags_sort


	def makeCluster(self, orient, chro, spec = None):
		if spec == None:
			spec = "human"

		res = []
		frags = self.sortFrags(orient, chro, spec)

		frag_clus = FragmentCluster.FragmentCluster()
		nfrags    = self.nfrags(orient, chro)

		for i in range(0, nfrags):
			# print("Fragment: %d of %d" %(i, nfrags))

			if frag_clus.isEmpty() or frag_clus.mergeStatus(frags[i], spec) == 1:
				frag_clus.addToCluster(frags[i])

				continue

			if frag_clus.mergeStatus(frags[i], spec) == 0:
				res.append(frag_clus)
				frag_clus = FragmentCluster.FragmentCluster()
				frag_clus.addToCluster(frags[i])

				continue

		res.append(frag_clus)
		return res


if __name__ == "__main__":
	import __ParseArgs
	args = __ParseArgs.parseArgs()
	args.infile = "tmp.sam"

	import pprint
	pp = pprint.PrettyPrinter(indent = 4)

	ideogram = IdeogramAlign(args)
	# pp.pprint(ideogram.data)

	pp.pprint(ideogram.getAllOrientation())
	pp.pprint(ideogram.getAllChromosome())
