#!/usr/bin/env python

# This package defines a class built by pairSamAlign object. The alignments of two reads with / without supplementary alignments were merged by figuring out their configuration, i.e., chromosome, strand and coordinates of left and right part of fragment. Therefore, the alignments of two reads were transformed into alignments of the fragment.

import re
import copy

import __Utils



##### Auxilary functions:
#

mergeSegs         = __Utils.mergeSegs

#
#####

"""
$ref = {
        rnam => "NA",
        left => {
            chro => "NA",
            strand => ".",
            start => 0,
            stop => 0,
            sequ => "NA",
            qual => "NA",
            rev => 0,
            circ => 0,
        },
        right => {
            chro => "NA",
            strand => ".",
            start => 0,
            stop => 0,
            sequ => "NA",
            qual => "NA",
            rev => 0,
            circ => 0,
        },
        orientation => "Human:+ | Virus:+",
        breakpoint_stat => 0,
        breakpoint => {
            human => {
                chro => "NA",
                strand => ".",
                start => 0,
                stop => 0,
            },
            virus => {
                chro => "NA",
                strand => ".",
                start => 0,
                stop => 0,
            },
        },
        homolog => {
            length => "NA",
            human => {
                chro => "NA",
                strand => "NA",
                start => 0,
                stop => 0,
            },
            virus => {
                chro => "NA",
                strand => "NA",
                start => 0,
                stop => 0,
            },
        },
		gap => {
			gap_flag   => 0,
			gap_chro   => "NA",
			gap_strand => ".",
			gap_start  => 0,
			gap_stop   => 0,
		},
        homo_discord => 0,     # set 1 if there are two Homolog records and
                               #they are discordant;
};
"""

class FragmentAlign:
    __slots__ = [
        "rnam", "left", "right",
        "breakpoint_stat", "homolog",
        "gap", "orientation", "breakpoint",
        "para",
    ]

    def __init__(self, ref, para):
        self.rnam            = ref["rnam"]
        self.left            = ref["left"]
        self.right           = ref["right"]
        self.breakpoint_stat = ref["breakpoint_stat"]
        self.homolog         = ref["homolog"]
        self.gap             = ref["gap"]
        
        # set 1 if there are two Homolog records and they are discordant;
        # self.homolog = ref["homolog"]


        if ref["rnam"] == "NA":
            return

        # del ref["stat"]

        chrVirus = para.virus

        strand_left  = ref["left"]["strand"]
        strand_right = ref["right"]["strand"]

        if ref["right"]["chro"] == chrVirus:
            self.orientation = "Human:%s | Virus:%s" %(strand_left, strand_right)
        else:
            self.orientation = "Virus:%s | Human:%s" %(strand_left, strand_right)

        breakpoint = dict(
            human = {
                "chro"   : "NA",
                "strand" : ".",
                "start"  : 0,
                "stop"   : 0,
            },
            virus = {
                "chro"   : "NA",
                "strand" : ".",
                "start"  : 0,
                "stop"   : 0,
            }
        )

        pat1 = re.compile("Human:+ | V")
        pat2 = re.compile("Virus:[+-] | H")

        if ref["breakpoint_stat"] == 0:
            if re.search(pat1, self.orientation) != None:
                breakpoint["human"]["chro"]   = ref["left"]["chro"]
                breakpoint["human"]["strand"] = ref["left"]["strand"]
                breakpoint["human"]["start"]  = ref["left"]["stop"]
                breakpoint["human"]["stop"]   = "NA"

                breakpoint["virus"]["chro"]   = chrVirus
                breakpoint["virus"]["strand"] = ref["right"]["strand"]
                breakpoint["virus"]["start"]  = "NA"
                breakpoint["virus"]["stop"]   = ref["right"]["start"]
            elif re.search(pat2, self.orientation) != None:
                breakpoint["human"]["chro"]   = ref["right"]["chro"]
                breakpoint["human"]["strand"] = ref["right"]["strand"]
                breakpoint["human"]["start"]  = "NA"
                breakpoint["human"]["stop"]   = ref["right"]["start"]

                breakpoint["virus"]["chro"]   = chrVirus
                breakpoint["virus"]["strand"] = ref["left"]["strand"]
                breakpoint["virus"]["start"]  = ref["left"]["stop"]
                breakpoint["virus"]["stop"]   = "NA"

            self.breakpoint = breakpoint
        elif ref["breakpoint_stat"] == 1:
            if re.search(pat1, self.orientation) != None:
                breakpoint["human"]["chro"]   = ref["left"]["chro"]
                breakpoint["human"]["strand"] = ref["left"]["strand"]
                breakpoint["human"]["start"]  = ref["left"]["stop"]
                breakpoint["human"]["stop"]   = ref["left"]["stop"]

                breakpoint["virus"]["chro"]   = chrVirus
                breakpoint["virus"]["strand"] = ref["right"]["strand"]
                breakpoint["virus"]["start"]  = ref["right"]["start"]
                breakpoint["virus"]["stop"]   = ref["right"]["start"]
            elif re.search(pat2, self.orientation) != None:
                breakpoint["human"]["chro"]   = ref["right"]["chro"]
                breakpoint["human"]["strand"] = ref["right"]["strand"]
                breakpoint["human"]["start"]  = ref["right"]["start"]
                breakpoint["human"]["stop"]   = ref["right"]["start"]

                breakpoint["virus"]["chro"]   = chrVirus
                breakpoint["virus"]["strand"] = ref["left"]["strand"]
                breakpoint["virus"]["start"]  = ref["left"]["stop"]
                breakpoint["virus"]["stop"]   = ref["left"]["stop"]

            self.breakpoint = breakpoint
        elif ref["breakpoint_stat"] == 2:
            self.breakpoint = ref["homolog"]
            # del self.breakpoint["length"]

        self.para = para


    def virusLength(self):
        return self.para.vlen


    def is_dummy(self):
        if self.rnam == "NA":
            return 1


    def chro(self, posi):
        # posi: left or right
        ss = "self.%s" %(posi)
        frag = eval(ss)

        return frag["chro"]


    def strand(self, posi):
        ss = "self.%s" %(posi)
        frag = eval(ss)

        return frag["strand"]


    def start(self, posi):
        ss = "self.%s" %(posi)
        frag = eval(ss)

        return frag["start"]


    def stop(self, posi):
        ss = "self.%s" %(posi)
        frag = eval(ss)

        return frag["stop"]


    def sequ(self, posi):
        ss = "self.%s" %(posi)
        frag = eval(ss)

        return frag["sequ"]


    def qual(self, posi):
        ss = "self.%s" %(posi)
        frag = eval(ss)

        return frag["qual"]


    def hasHomolog(self):
        if self.homolog["length"] > 0:
            return 1
        else:
            return 0


    def revStat(self, posi):
        ss = "self.%s" %(posi)
        frag = eval(ss)

        return frag["rev"]


    def circStat(self, posi):
        ss = "self.%s" %(posi)
        frag = eval(ss)

        return frag["circ"]


    # def orientation(self):
    #     return self.orientation


    # def breakpoint(self):
    #     return self.breakpoint


    # def breakpoint_stat(self):
    #     return self.breakpoint_stat


    def configPosi(self):
        res = dict(
            human = "NA",
            virus = "NA",
        )

        chrVirus = self.para.virus

        if self.chro("left") == chrVirus:
            res["human"] = "right"
            res["virus"] = "left"
        else:
            res["human"] = "left"
            res["virus"] = "right"
        
        return res


    def humanPart(self):
        posi = self.configPosi()

        ss = "self.%s" %(posi["human"])
        res = eval(ss)
        return res


    def virusPart(self):
        posi = self.configPosi()

        ss = "self.%s" %(posi["virus"])
        res = eval(ss)
        return res


    def getCase(self):
        chr_left  = self.chro("left")
        chr_right = self.chro("right")
        cha_left  = self.strand("left")
        cha_right = self.strand("right")

        return "%s:%s | %s:%s" %(chr_left, cha_left, chr_right, cha_right)


    def sameConfig(self, obje):
        case1 = self.getCase()
        case2 = obje.getCase()

        if case1 == case2:
            return 1
        else:
            return 0


    def sameAlignCord(frag1, frag2):
        #  1: same
        #  0: different
        # -1: with different configuration

        if not frag1.sameConfig(frag2):
            return -1

        sta1_left  = frag1.start("left")
        sto1_left  = frag1.stop("left")
        sta1_right = frag1.start("right")
        sto1_right = frag1.stop("right")

        sta2_left  = frag2.start("left")
        sto2_left  = frag2.stop("left")
        sta2_right = frag2.start("right")
        sto2_right = frag2.stop("right")

        sta1_gap = frag1.gap["gap_start"]
        sto1_gap = frag1.gap["gap_stop"]
        sta2_gap = frag2.gap["gap_start"]
        sto2_gap = frag2.gap["gap_stop"]

        # Homolog
        sta1_homo_human = frag1.homolog["human"]["start"]
        sto1_homo_human = frag1.homolog["human"]["stop"]
        sta2_homo_human = frag2.homolog["human"]["start"]
        sto2_homo_human = frag2.homolog["human"]["stop"]

        sta1_homo_virus = frag1.homolog["virus"]["start"]
        sto1_homo_virus = frag1.homolog["virus"]["stop"]
        sta2_homo_virus = frag2.homolog["virus"]["start"]
        sto2_homo_virus = frag2.homolog["virus"]["stop"]

        if sta1_left == sta2_left and sto1_left == sto2_left and sta1_right == sta2_right and sto1_right == sto2_right:
            if sta1_gap == sta2_gap and sto1_gap == sto2_gap and sta1_homo_human == sta2_homo_human and sto1_homo_human == sto2_homo_human and sta1_homo_virus == sta2_homo_virus and sto1_homo_virus == sto2_homo_virus:
                return 1
            else:
                return 0
        else:
            return 0


    def fragStr(self):
        sta_left  = self.start("left")
        sto_left  = self.stop("left")
        sta_right = self.start("right")
        sto_right = self.stop("right")

        sta_gap = self.gap["gap_start"]
        sto_gap = self.gap["gap_stop"]

        # Homolog
        sta_homo_human = self.homolog["human"]["start"]
        sto_homo_human = self.homolog["human"]["stop"]
        sta_homo_virus = self.homolog["virus"]["start"]
        sto_homo_virus = self.homolog["virus"]["stop"]

        sss = f"{sta_left} - {sto_left}; {sta_right} - {sto_right}; {sta_gap} - {sto_gap}; {sta_homo_human} - {sto_homo_human}; {sta_homo_virus} - {sto_homo_virus}"
        return sss


    def combineFrags(frag1, frag2):
        res = copy.deepcopy(frag2)

        # import pprint
        # pp = pprint.PrettyPrinter(indent = 4)

        # print("Frag 1: " + frag1.rnam)
        # pp.pprint(frag1.left)
        # pp.pprint(frag1.right)

        # print("Frag 2: " + frag2.rnam)
        # pp.pprint(frag2.left)
        # pp.pprint(frag2.right)

        seq1_left  = frag1.sequ("left")
        qua1_left  = frag1.qual("left")
        seq1_right = frag1.sequ("right")
        qua1_right = frag1.qual("right")

        # print("Seq 1 Left: \n" + seq1_left)
        # print(qua1_left)

        seq2_left  = frag2.sequ("left")
        qua2_left  = frag2.qual("left")
        seq2_right = frag2.sequ("right")
        qua2_right = frag2.qual("right")

        # print("Seq 2 Left: \n" + seq2_left)
        # print(qua2_left)


        # Left sequence merge:
        nleft = len(seq1_left)
        sequ_left = ""
        qual_left = ""
        for i in range(0, nleft):
            # print("Cord: " + str(i))

            bas1 = seq1_left[i]
            qua1 = qua1_left[i]
            bas2 = seq2_left[i]
            qua2 = qua2_left[i]

            if bas1 == "-":
                sequ_left += "-"
                qual_left += "-"
                continue
            
            if qua1 >= qua2:
                sequ_left += bas1
                qual_left += qua1
                continue
            else:
                sequ_left += bas2
                qual_left += qua2
                continue
        
        res.left["sequ"] = sequ_left
        res.left["qual"] = qual_left
        
        # Right sequence merge:
        nright = len(seq1_right)
        sequ_right = ""
        qual_right = ""
        for i in range(0, nright):
                bas1 = seq1_right[i]
                qua1 = qua1_right[i]
                bas2 = seq2_right[i]
                qua2 = qua2_right[i]
                
                if bas1 == "-":
                    sequ_right += "-"
                    qual_right += "-"
                    continue
                
                if qua1 >= qua2:
                    sequ_right += bas1
                    qual_right += qua1
                    continue
                else:
                    sequ_right += bas2
                    qual_right += qua2
                    continue
        
        res.right["sequ"] = sequ_right
        res.right["qual"] = qual_right
        
        return res
    

    def figurePosi(frag1, frag2, spec = "human"):
        if frag1.sameConfig(frag2) == 0:
            return "NA", "NA"
        
        config = frag1.configPosi()
        sta1 = frag1.start( config(spec) )
        sto1 = frag1.stop( config(spec) )
        sta2 = frag2.start( config(spec) )
        sto2 = frag2.stop( config(spec) )

        if frag1[ config(spec) ]["rev"] == 0:
            if sta1 < sta2:
                return frag1, frag2
            elif sta1 == sta2:
                if sto1 <= sto2:
                    return frag1, frag2
                else:
                    return frag2, frag1
            else:
                return frag2, frag1
        
        if frag1[ config(spec) ]["rev"] == 1:
            if sta1 > sta2:
                return frag1, frag2
            elif sta1 == sta2:
                if sto1 >= sto2:
                    return frag1, frag2
                else:
                    return frag2, frag1
            else:
                return frag2, frag1
    

    # Merge fragments, designed to merge the sequences of two fragments
    def mergeFrags(frag1, frag2, circular_stat):
        vlen = frag1.virusLength()
        res = copy.deepcopy(frag1)

        fragL, gapL = mergeSegs(frag1.left, frag2.left, vlen, circular_stat)
        res.left = fragL

        fragR, gapR = mergeSegs(frag1.right, frag2.right, vlen, circular_stat)
        res.right = fragR

        return res
    

    def figureBreakpoint(self, spec = "human"):
        config = self.configPosi()
        posi   = config[ spec ]
        sta    = self.breakpoint[spec]["start"]
        sto    = self.breakpoint[spec]["stop"]

        res = dict(
            chro   = self.chro(posi),
            strand = self.strand(posi),
            cord   = 0,
            sta    = sta,
            sto    = sto,
        )

        if self.breakpoint_stat == 0:
            if sto == "NA":
                res["cord"] = sta
            
            if sta == "NA":
                res["cord"] = sto
        else:
            res["cord"] = (sta + sto)/2
        
        return res
    

    def strFragmentAlign(self):
        chrVirus = self.para.virus

        rnam    = self.rnam
        chrL    = self.chro("left")
        strandL = self.strand("left")
        staL    = self.start("left")
        stoL    = self.stop("left")

        chrR    = self.chro("right")
        strandR = self.strand("right")
        staR    = self.start("right")
        stoR    = self.stop("right")

        breakpoint = self.breakpoint

        bp_tmp = None
        if self.chro("left") == chrVirus:
            bp_left = "%s:%s:%s:%s" %(
                breakpoint["virus"]["chro"], 
                breakpoint["virus"]["strand"], 
                str(breakpoint["virus"]["start"]), 
                str(breakpoint["virus"]["stop"])
            )
            bp_right = "%s:%s:%s:%s" %(
                breakpoint["human"]["chro"],
                breakpoint["human"]["strand"],
                str(breakpoint["human"]["start"]),
                str(breakpoint["human"]["stop"])
            )
            bp_tmp = "%s | %s" %(bp_left, bp_right)
        else:
            bp_left = "%s:%s:%s:%s" %(
                breakpoint["human"]["chro"],
                breakpoint["human"]["strand"],
                str(breakpoint["human"]["start"]),
                str(breakpoint["human"]["stop"])
            )
            bp_right = "%s:%s:%s:%s" %(
                breakpoint["virus"]["chro"],
                breakpoint["virus"]["strand"],
                str(breakpoint["virus"]["start"]),
                str(breakpoint["virus"]["stop"])
            )
            bp_tmp = "%s | %s" %(bp_left, bp_right)
        
        bp_type = None
        if self.breakpoint_stat == 0:
            bp_type = 0
        else:
            bp_type = 1
        
        vec = [
            rnam, chrL, strandL, str(staL), str(stoL), chrR,
            strandR, str(staR), str(stoR), bp_tmp, str(bp_type)
        ]
        
        return "\t".join(vec)
    

    """
    Returned values:
        -1: circular;
         0: head;
         1: tail;
    """
    def headOrTail(self):
        posiV = self.configPosi()["virus"]

        if self.circStat(posiV) == 1:
            return -1
        
        sta = self.start(posiV)
        sto = self.stop(posiV)
        cod = (sta + sto)/2

        if self.para.vlen - cod < cod - 1:
            return 1
        else:
            return 0
    

    def virusFragCirc(self, frag):
        flag1 = self.headOrTail()
        flag2 = frag.headOrTail()

        if flag1 ^ flag2:
            return 1
        else:
            return 0
    

    # Test if fragment of virus is too long (>= 500 bp), which is abnormal
    def virusFragLengthTest(self):
        posiV = self.configPosi("virus")

        sta = self.start(posiV)
        sto = self.stop(posiV)

        dis = None
        if self.circStat(posiV) == 0:
            dis = abs(sta - sto)
        else:
            if self.revStat(posiV) == 0:
                dis = (self.virusLength() - sta + 1) + sto
            else:
                dis = (self.virusLength() - sto + 1) + sta
        
        if dis >= 500:
            return 1
        else:
            return 0
    

if __name__ == "__main__":
    import __ParseArgs
    import SingleSamAlign
    import SingleAlign
    import PairSamAlign

    import pprint
    pp = pprint.PrettyPrinter(indent = 4)

    args = __ParseArgs.parseArgs()

    rec = {}
    reads = []

    with open("Test_Data_SingleAlign.sam", "r") as f:
        for line in f:
            tmp = line.split("\t")

            idx = tmp[0].find("_")
            nam = tmp[0][:idx]

            if not nam in rec:
                reads.append(nam)
                rec[ nam ] = []

            rec[ nam ].append(line[:-1])

    for read_tmp in reads:
        id = read_tmp
        print(">%s" %(id))

        pairAln = PairSamAlign.PairSamAlign(rec[read_tmp], args)
        frag_ref = pairAln.toFragment()
        frag = FragmentAlign(frag_ref, args)

        print(isinstance(frag, FragmentAlign))
        print(frag.fragStr())

        # pp.pprint(frag)
        pp.pprint(frag_ref)

        frag_tmp = copy.deepcopy(frag)
        frag.breakpoint_stat = 5
        print("Breakpoint status for frag_tmp: %d" %(frag_tmp.breakpoint_stat))

        config = frag.configPosi()
        pp.pprint(config)

        ss = frag.getCase()
        print(ss)

        bp_human = frag.figureBreakpoint("human")
        pp.pprint(bp_human)
        
        bp_virus = frag.figureBreakpoint("virus")
        pp.pprint(bp_virus)

        strFrag = frag.strFragmentAlign()
        print(strFrag)

        flag = frag.headOrTail()
        print("Head or Tail: %d" %(flag))

        del frag.rnam
        try:
            print("read name: ", frag.rnam)
        except AttributeError:
            print("read name has been deleted from frag.")
        finally:
            print("Excuting final step.")

        print("\n ---------- <END> ----------\n")

        break

