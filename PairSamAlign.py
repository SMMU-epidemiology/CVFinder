#!/usr/bin/env python

import re
import pprint
import copy

import __Utils
import SingleSamAlign
import SingleAlign

##### Auxilary functions:
#

compSeq           = __Utils.compSeq
pruneSequ         = __Utils.pruneSequ
revCompSequ       = __Utils.revCompSequ
set_reverse_cigar = __Utils.set_reverse_cigar
mergeSegs         = __Utils.mergeSegs

# bioExplainCigar = __Utils.bioExplainCigar

#
#####

pp = pprint.PrettyPrinter(indent = 4)

class PairSamAlign:
    __slots__ = [
        "read1", "read1_sa", "read1_circ",
        "read2", "read2_sa", "read2_circ",
        "para", "sign",
    ]
    
    def __init__(self, ref, para):
        self.read1      = ""
        self.read1_sa   = ""
        self.read1_circ = 0
        self.read2      = ""
        self.read2_sa   = ""
        self.read2_circ = 0
        self.para       = para
        self.sign       = 0

        for ss in ref:
            read_tmp = SingleSamAlign.SingleSamAlign(ss, para)

            if read_tmp.is_OK():
                id = read_tmp.idName()

                ss = "self.%s" %(id)
                aln = eval(ss)
                if aln == "":
                    if id == "read1":
                        self.read1 = read_tmp

                    if id == "read1_sa":
                        self.read1_sa = read_tmp

                    if id == "read2":
                        self.read2 = read_tmp

                    if id == "read2_sa":
                        self.read2_sa = read_tmp
                else:
                    self.sign = 1

        if self.read1 == "":
            self.read1 = SingleSamAlign.SingleSamAlign("", para)

        if self.read1_sa == "":
            self.read1_sa = SingleSamAlign.SingleSamAlign("", para)

        if self.read2 == "":
            self.read2 = SingleSamAlign.SingleSamAlign("", para)

        if self.read2_sa == "":
            self.read2_sa = SingleSamAlign.SingleSamAlign("", para)

        #
        if self.sign == 1:
            return


        if not self.read1.is_dummy() and self.read1.is_not_simple():
            self.read1.fixCigar()

        if not self.read2.is_dummy() and self.read2.is_not_simple():
            self.read2.fixCigar()

        if not self.read1_sa.is_dummy():
            read1 = self.read1
            read1_sa = self.read1_sa

            if read1_sa.is_not_simple():
                read1_sa.fixCigar()

            if not read1.args.linear:
                circ = read1.tailToHead(read1_sa)

                if circ["stat"] == 1:
                    self.read1 = circ["merge"]
                    self.read1_sa = SingleSamAlign.SingleSamAlign("", para)
                    self.read1_circ = 1

        if not self.read2_sa.is_dummy():
            read2 = self.read2
            read2_sa = self.read2_sa

            if read2_sa.is_not_simple():
                read2_sa.fixCigar()

            if not read2.args.linear:
                circ = read2.tailToHead(read2_sa)

                if circ["stat"] == 1:
                    self.read2 = circ["merge"]
                    self.read2_sa = SingleSamAlign.SingleSamAlign("", para)
                    self.read2_circ = 1


    def vlen(self):
        return self.para.vlen


    def rnam(self):
        return self.read1.rnam


    # "0" is OK
    def mapqCheck(self, thre = 30):
        if self.read1.mapq < thre or self.read1_sa.mapq < thre or self.read2.mapq < thre or self.read2_sa.mapq < thre:
            return 1
        else:
            return 0


    # "0" is OK
    def flagCheck(self):
        if self.sign == 0:
            return 0
        else:
            return 1


    # "0" is OK
    def is_unmapped(self):
        if self.read1.flag & 4 or self.read1_sa.flag & 4 or self.read2.flag & 4 or self.read2_sa.flag & 4:
            return 1
        else:
            return 0


    # "0" is OK
    def is_discordant(self):
        chr1 = self.read1.chro
        chr1_sa = self.read1_sa.chro
        chr2 = self.read2.chro
        chr2_sa = self.read2_sa.chro

        ref = {}
        if chr1 != "NA":
            ref[chr1] = None

        if chr1_sa != "NA":
            ref[chr1_sa] = None

        if chr2 != "NA":
            ref[chr2] = None

        if chr2_sa != "NA":
            ref[chr2_sa] = None

        keys = ref.keys()
        if len(keys) > 2 or len(keys) == 1:
            return 1
        else:
            return 0


    # Test if the alignment is pair mapped;
    # 0: OK; 1: not OK.
    def not_pair_mapped(self):
        n = 0

        if self.read1.chro != "NA":
            n += 1

        if self.read1_sa.chro != "NA":
            n += 1

        if self.read2.chro != "NA":
            n += 1

        if self.read2_sa.chro != "NA":
            n += 1

        if n < 2:
            return 1
        else:
            return 0


    def is_dummy(self):
        if self.read1_sa.is_dummy() and self.read2_sa.is_dummy():
            return 1
        else:
            return 0


    def singleEndSeq(self):
        if self.read2.is_dummy() and self.read2_sa.is_dummy():
            if self.read1.isSingleEndSeq() and self.read1_sa.isSingleEndSeq():
                return 1
            else:
                return 0
        else:
            return 0


    def is_strand_abnormal(self):
        if self.is_dummy():
            return 0

        res = {
            "Human" : [],
            "Virus" : [],
        }
        spec = None

        spec = self.read1.getSpec()
        res[spec].append(self.read1.strand())

        if not self.read1_sa.is_dummy():
            spec = self.read1_sa.getSpec()
            res[spec].append(self.read1_sa.strand())

        spec = self.read2.getSpec()

        # KeyError Bug Fixed:
        if not spec in res:
            res[spec] = []

        res[spec].append(self.read2.strand())

        if not self.read2_sa.is_dummy():
            spec = self.read2_sa.getSpec()
            res[spec].append(self.read2_sa.strand())

        flagH = None
        flagV = None
        if len(res["Human"]) == 1:
            flagH = 0
        else:
            if res["Human"][0] == res["Human"][1]:
                flagH = 1
            else:
                flagH = 0

        if len(res["Virus"]) == 1:
            flagV = 0
        else:
            if res["Virus"][0] == res["Virus"][1]:
                flagV = 1
            else:
                flagV = 0

        if flagH == 1 or flagV == 1:
            return 1
        else:
            return 0


    def alignQC(self, thre = 30):
        if self.flagCheck() == 1:
            return 1

        if self.mapqCheck(thre):
            return 1

        if self.is_unmapped():
            return 1

        if self.is_discordant():
            return 1

        if self.not_pair_mapped():
            return 1

        if self.is_strand_abnormal():
            return 1

        return 0


    def config_OK(self):
        chrVirus = self.para.virus
        para = self.read1.args

        # print("\n -----  Parameters: -----")
        # print(para)

        # if para.linear == False:
        #   print("The virus is circular.")

        # print(" -----  /////////// -----\n")
        
        """
        my $res = {
            config_stat: 1: OK; 0: Not OK
            act: rev: 0 / 1; comp: 0 / 1; strand_comp: 0 / 1;
            left: SingleAlignRef / config;
            right: SingleAlignRef / config;
            left_type: SingleAlignRef / config;
            right_type: SingleAlignRef / config;
            left_id: read1 / read2;
            right_id: read1 / read2;
        };

        """

        res = {
            # 1: OK; 0: Not OK  
            "config_stat" : 0,  
            "act"         : "NA",
            "left"        : "NA",
            "right"       : "NA",
            "left_type"   : "NA",
            "right_type"  : "NA",
            "left_id"     : "NA",
            "right_id"    : "NA",
        }

        # ==================================================
        # The virus "Reverse" operation in SingleAlign objects has been done by alignConfiguration function in SingleAlign.py during the process of cigar reverse. Therefore some HBV "rev" flags are always "0", which isn't a bug!
        # ==================================================

        act = {
            "read1" : {
                "rev"         : 0,
                "comp"        : 0,
                "strand_comp" : 0,
            },
            "read1_sa" : {
                "rev"         : 0,
                "comp"        : 0,
                "strand_comp" : 0,
            },
            "read2" : {
                "rev"         : 0,
                "comp"        : 0,
                "strand_comp" : 0,
            },
            "read2_sa" : {
                "rev"         : 0,
                "comp"        : 0,
                "strand_comp" : 0,
            },
        }

        res["act"] = act

        if self.is_dummy():
            res["config_stat"] = 1
            return res

        flag1 = self.read1_sa.is_dummy()
        flag2 = self.read2_sa.is_dummy()

        ref = {
            "primary" : "",
            "supple"  : "",
        }
        config    = None
        human     = None
        virus     = None
        virus_rev = None
        this      = None
        other     = None

        ###################################################
        ##             Single End Sequencing
        ###################################################

        if self.singleEndSeq():
            ref["primary"] = self.read1
            ref["supple"]  = self.read1_sa
            this           = "read1"

            align = SingleAlign.SingleAlign(ref)

            # "is_hs_virus" returns 1 if there's a problem.
            if align.is_hs_virus() == 1:
                return res

            config = align.alignConfiguration()

            if config["left"]["chro"] == chrVirus:
                human     = "right"
                virus     = "left"
                virus_rev = "left_rev"
            else:
                human     = "left"
                virus     = "right"
                virus_rev = "right_rev"

            res["config_stat"] = 1

            if config[human]["strand"] == "+":
                res["left"]      = config
                res["left_type"] = "config"
                res["left_id"]   = this

                if human == "left":
                    # Human: Left +;
                    if config[virus_rev] == 0:
                        # Human: Left +; Virus: no reverse

                        return res
                    else:
                        # Human: Left +; Virus: reverse

                        virus_id = this
                        if config["right_id"] == "supple":
                            virus_id += "_sa"
                            res["act"][virus_id]["comp"] = 1

                        return res
                else:
                    # Human: Right +; Virus: NA

                    virus_id = this
                    if config[virus_rev] == 1:
                        if config["left_id"] == "supple":
                            virus_id += "_sa"

                        res["act"][virus_id]["comp"] = 1

                    return res
            else:
                res["left"]      = config
                res["left_type"] = "config"
                res["left_id"]   = this

                if human == "right":
                    # Human: Right -;

                    human_id = None
                    virus_id = None
                    if config["right_id"] == "supple":
                        human_id = "%s_sa" %(this)
                        virus_id = this
                    else:
                        human_id = this
                        virus_id = "%s_sa" %(this)

                    if config[virus_rev] == 0:
                        # Human: Right -; Virus: no reverse

                        res["act"][human_id]["strand_comp"] = 1
                        res["act"][virus_id]["strand_comp"] = 1

                        return res
                    else:
                        # Human: Right -; Virus: reverse

                        res["act"][human_id]["strand_comp"] = 1
                        res["act"][virus_id]["comp"]        = 1
                        res["act"][virus_id]["strand_comp"] = 1

                        return res
                else:
                    # Human: Left -; Virus: NA

                    human_id = None
                    virus_id = None
                    if config["left_id"] == "supple":
                        human_id = "%s_sa" %(this)
                        virus_id = this
                    else:
                        human_id = this
                        virus_id = "%s_sa" %(this)

                    res["act"][human_id]["strand_comp"] = 1
                    res["act"][virus_id]["strand_comp"] = 1

                    if config[virus_rev] == 1:
                        res["act"][virus_id]["comp"] = 1

                    return res


        # import pprint
        # pp = pprint.PrettyPrinter(indent = 4)

        ###################################################
        ##             Pair End Sequencing
        ###################################################

        if flag1 == 0 and flag2 == 1:
            ref["primary"] = self.read1
            ref["supple"]  = self.read1_sa
            this           = "read1"
            other          = "read2"
        elif flag1 == 1 and flag2 == 0:
            ref["primary"] = self.read2
            ref["supple"]  = self.read2_sa
            this           = "read2"
            other          = "read1"

        if flag1 != 0 or flag2 != 0:
            # print("----- Line: 494 -----\n")

            align = SingleAlign.SingleAlign(ref)

            # print(">SingleAlign Construction:")
            # print(">Primary Alignment:")
            # pp.pprint(align.primaryAlignment())
            # print(">Supplementary Alignment:")
            # pp.pprint(align.suppleAlignment())


            # "is_hs_virus" returns 1 if there's a problem.
            if align.is_hs_virus():
                return res

            config = align.alignConfiguration()

            # print("\n>alignConfiguration:")
            # pp.pprint(config)


            if config["left"]["chro"] == "chrVirus":
                human     = "right"
                virus     = "left"
                virus_rev = "left_rev"
            else:
                human     = "left"
                virus     = "right"
                virus_rev = "right_rev"

            if config[human]["strand"] == "+":
                res["left"]       = config
                res["left_type"]  = "config"
                res["right"]      = eval("self.%s.toSingleAlignRef()" %(other))
                res["right_type"] = "SingleAlignRef"

                res["left_id"]    = this
                res["right_id"]   = other

                if human == "left":
                    # Human: Left +;

                    # print("----- Line: 525 -----\n")

                    ss = "self.%s" %(other)
                    aln_tmp = eval(ss)

                    if aln_tmp.chro != chrVirus:
                        return res
                    
                    if config[virus_rev] == 0:
                        # Human: Left +; Virus: no reverse;

                        # print("----- Line: 536 -----\n")

                        sta_other = aln_tmp.toSingleAlignRef()["start"]
                        sto_other = aln_tmp.toSingleAlignRef()["stop"]
                        cod_other = (abs(sta_other) + abs(sto_other))/2

                        sta_suppl = config[virus]["start"]
                        sto_suppl = config[virus]["stop"]
                        cod_suppl = (sta_suppl + sto_suppl)/2

                        if aln_tmp.toSingleAlignRef()["circ"] == 1:
                            # print("----- Line: 541 -----\n")

                            if abs(sta_other) >= sta_suppl:
                                res["config_stat"] = 1
                                res["act"][other]["strand_comp"] = 1
                        else:
                            dis_linear = abs(cod_other - cod_suppl)
                            dis_circular = self.vlen() - cod_suppl + 1 + cod_other

                            if dis_linear <= dis_circular:
                                # print("----- Line: 551 -----\n")

                                if sta_other >= sta_suppl:
                                    res["config_stat"] = 1
                                    res["act"][other]["strand_comp"] = 1
                            else:
                                # print("----- Line: 557 -----\n")

                                if cod_other < cod_suppl:
                                    res["config_stat"] = 1
                                    res["act"][other]["strand_comp"] = 1
                        
                        return res
                    else:
                        # Human: Left +; Virus: reverse;

                        # print("----- Line: 573 -----\n")

                        ss = "self.%s" %(other)
                        aln_tmp = eval(ss)

                        sta_other = aln_tmp.toSingleAlignRef()["stop"]
                        sto_other = aln_tmp.toSingleAlignRef()["start"]
                        cod_other = (abs(sta_other) + abs(sto_other))/2

                        sta_suppl = config[virus]["start"]
                        sto_suppl = config[virus]["stop"]
                        cod_suppl = (sta_suppl + sto_suppl)/2

                        if aln_tmp.toSingleAlignRef()["circ"] == 1:
                            # print("----- Line: 579 -----\n")

                            if abs(sta_other) <= sta_suppl:
                                res["config_stat"] = 1

                                virus_id = this
                                if config["right_id"] == "supple":
                                    virus_id += "_sa"

                                res["act"][virus_id]["comp"] = 1
                                res["act"][other]["rev"] = 1
                                res["act"][other]["comp"] = 1
                                res["act"][other]["strand_comp"] = 1
                        else:
                            if para.linear == True:
                                # print("----- Line: 594 -----\n")

                                if sta_other <= sta_suppl:
                                    res["config_stat"] = 1

                                    virus_id = this
                                    if config["right_id"] == "supple":
                                        virus_id += "_sa"

                                    res["act"][virus_id]["comp"] = 1
                                    res["act"][other]["rev"] = 1
                                    res["act"][other]["comp"] = 1
                                    res["act"][other]["strand_comp"] = 1
                            else:
                                dis_linear = abs(cod_other - cod_suppl)
                                dis_circular = self.vlen() - cod_other + 1 + cod_suppl

                                if dis_linear <= dis_circular:
                                    # print("----- Line: 612 -----\n")

                                    if sta_other <= sta_suppl:
                                        res["config_stat"] = 1

                                        virus_id = this
                                        if config["right_id"] == "supple":
                                            virus_id += "_sa"

                                        res["act"][virus_id]["comp"] = 1
                                        res["act"][other]["rev"] = 1
                                        res["act"][other]["comp"] = 1
                                        res["act"][other]["strand_comp"] = 1
                                else:
                                    # print("----- Line: 626 -----\n")

                                    if cod_other > cod_suppl:
                                        res["config_stat"] = 1

                                        virus_id = this
                                        if config["right_id"] == "supple":
                                            virus_id += "_sa"

                                        res["act"][virus_id]["comp"] = 1
                                        res["act"][other]["rev"] = 1
                                        res["act"][other]["comp"] = 1
                                        res["act"][other]["strand_comp"] = 1

                        return res
                else:
                    # Human: Right +; Virus: NA

                    ss = "self.%s" %(other)
                    aln_tmp = eval(ss)

                    if aln_tmp.chro == chrVirus:
                        return res

                    sta_other = aln_tmp.cord
                    sta_suppl = config[human]["start"]
                    sto_suppl = config[human]["stop"]

                    if abs(sta_other - sta_suppl) >= 2000:
                        # print("----- Line: 655 -----\n")

                        res["config_stat"] = 0
                        return res

                    # print("----- Line: 660 -----\n")
                    # print("sta_other: %d" %(sta_other))
                    # print("sta_suppl: %d" %(sta_suppl))

                    if sta_other >= sta_suppl:
                        # print("----- Line: 665 -----\n")

                        res["config_stat"] = 1
                        res["act"][other]["strand_comp"] = 1

                    virus_id = this
                    if config[virus_rev] == 1:
                        if config["left_id"] == "supple":
                            virus_id += "_sa"
                        res["act"][virus_id]["comp"] = 1

                    return res
            else:
                ss = "self.%s" %(other)
                aln_tmp = eval(ss)

                res["left"]       = aln_tmp.toSingleAlignRef()
                res["left_type"]  = "SingleAlignRef"
                res["right"]      = config
                res["right_type"] = "config"

                res["left_id"]  = other
                res["right_id"] = this

                if human == "right":
                    # Human: Right -;

                    if aln_tmp.chro != chrVirus:
                        return res

                    human_id = None
                    virus_id = None
                    if config["right_id"] == "supple":
                        human_id = "%s_sa" %(this)
                        virus_id = this
                    else:
                        human_id = this
                        virus_id = "%s_sa" %(this)

                    if config[virus_rev] == 0:
                        # Human: Right -; Virus: no reverse;

                        sta_other = aln_tmp.toSingleAlignRef()["start"]
                        sto_other = aln_tmp.toSingleAlignRef()["stop"]
                        cod_other = (abs(sta_other) + abs(sto_other))/2

                        sta_suppl = config[virus]["start"]
                        sto_suppl = config[virus]["stop"]
                        cod_suppl = (sta_suppl + sto_suppl)/2

                        if aln_tmp.toSingleAlignRef()["circ"] == 1:
                            # print("----- Line: 712 -----\n")

                            if abs(sto_other) <= sto_suppl:
                                res["config_stat"] = 1
                                res["act"][human_id]["strand_comp"] = 1
                                res["act"][virus_id]["strand_comp"] = 1
                        else:
                            if para.linear == True:
                                # print("----- Line: 720 -----\n")

                                if abs(sto_other) <= sto_suppl:
                                    res["config_stat"] = 1
                                    res["act"][human_id]["strand_comp"] = 1
                                    res["act"][virus_id]["strand_comp"] = 1
                            else:
                                dis_linear = abs(cod_other - cod_suppl)
                                dis_circular = self.vlen() - cod_other + 1 + cod_suppl

                                if dis_linear <= dis_circular:
                                    # print("----- Line: 731 -----\n")

                                    if abs(sto_other) <= sto_suppl:
                                        res["config_stat"] = 1
                                        res["act"][human_id]["strand_comp"] = 1
                                        res["act"][virus_id]["strand_comp"] = 1
                                else:

                                    # print("----- Line: 739 -----\n")

                                    if cod_other > cod_suppl:
                                        res["config_stat"] = 1
                                        res["act"][human_id]["strand_comp"] = 1
                                        res["act"][virus_id]["strand_comp"] = 1

                        return res
                    else:
                        # Human: Right -; Virus: reverse

                        sta_other = aln_tmp.toSingleAlignRef()["stop"]
                        sto_other = aln_tmp.toSingleAlignRef()["start"]
                        cod_other = (abs(sta_other) + abs(sto_other))/2

                        sta_suppl = config[virus]["start"]
                        sto_suppl = config[virus]["stop"]
                        cod_suppl = (sta_suppl + sto_suppl)/2

                        if aln_tmp.toSingleAlignRef()["circ"] == 1:
                            # print("----- Line: 759 -----\n")

                            if abs(sto_other) >= sto_suppl:
                                res["config_stat"] = 1

                                res["act"][human_id]["strand_comp"] = 1
                                res["act"][virus_id]["comp"] = 1
                                res["act"][virus_id]["strand_comp"] = 1
                                res["act"][other]["rev"] = 1
                                res["act"][other]["comp"] = 1
                        else:
                            if para.linear == True:
                                # print("----- Line: 771 -----\n")

                                if abs(sto_other) >= sto_suppl:
                                    res["config_stat"] = 1
                                    
                                    res["act"][human_id]["strand_comp"] = 1
                                    res["act"][virus_id]["comp"] = 1
                                    res["act"][virus_id]["strand_comp"] = 1
                                    res["act"][other]["rev"] = 1
                                    res["act"][other]["comp"] = 1
                            else:
                                dis_linear = abs(cod_other - cod_suppl)
                                dis_circular = self.vlen() - cod_suppl + 1 + cod_other

                                if dis_linear <= dis_circular:
                                    # print("----- Line: 786 -----\n")

                                    if abs(sto_other) >= sto_suppl:
                                        res["config_stat"] = 1

                                        res["act"][human_id]["strand_comp"] = 1
                                        res["act"][virus_id]["comp"] = 1
                                        res["act"][virus_id]["strand_comp"] = 1
                                        res["act"][other]["rev"] = 1
                                        res["act"][other]["comp"] = 1
                                else:
                                    # print("----- Line: 797 -----\n")

                                    if cod_other < cod_suppl:
                                        res["config_stat"] = 1

                                        res["act"][human_id]["strand_comp"] = 1
                                        res["act"][virus_id]["comp"] = 1
                                        res["act"][virus_id]["strand_comp"] = 1
                                        res["act"][other]["rev"] = 1
                                        res["act"][other]["comp"] = 1

                        return res
                else:
                    # Human: Left -; Virus: NA

                    ss = "self.%s" %(other)
                    aln_tmp = eval(ss)

                    if aln_tmp.chro == chrVirus:
                        return res

                    sto_other = aln_tmp.toSingleAlignRef()["stop"]
                    sta_suppl = config[human]["start"]
                    sto_suppl = config[human]["stop"]

                    if abs(sto_other - sto_suppl) >= 2000:
                        # print("----- Line: 823 -----\n")

                        res["config_stat"] = 0
                        return res

                    human_id = None
                    virus_id = None
                    if config["left_id"] == "supple":
                        human_id = "%s_sa" %(this)
                        virus_id = this
                    else:
                        human_id = this
                        virus_id = "%s_sa" %(this)

                    # print("----- Line: 837 -----\n")

                    if sto_other <= sto_suppl:
                        res["config_stat"] = 1
                        res["act"][human_id]["strand_comp"] = 1
                        res["act"][virus_id]["strand_comp"] = 1

                    if config[virus_rev] == 1:
                        res["act"][virus_id]["comp"] = 1

                    return res

        if flag1 == 0 and flag2 == 0:
            ref["primary"] = self.read1
            ref["supple"]  = self.read1_sa

            align1 = SingleAlign.SingleAlign(ref)

            if align1.is_hs_virus():
                return res

            config1 = align1.alignConfiguration()

            ref["primary"] = self.read2
            ref["supple"]  = self.read2_sa

            align2 = SingleAlign.SingleAlign(ref)

            if align2.is_hs_virus():
                return res

            config2 = align2.alignConfiguration()

            if config1["left"]["chro"] != config2["left"]["chro"] or config1["right"]["chro"] != config2["right"]["chro"]:
                return res

            if config1["left"]["stop"] != config2["left"]["stop"] or config1["right"]["start"] != config2["right"]["start"]:
                return res

            # print("----- Line: 876 -----\n")

            res["config_stat"] = 1
            res["left_type"]   = "config"
            res["right_type"]  = "config"

            if align1.human_strand() == "+":
                res["left"]     = config1
                res["left_id"]  = "read1"
                res["right"]    = config2
                res["right_id"] = "read2"

                res["act"]["read2"]["strand_comp"]    = 1
                res["act"]["read2_sa"]["strand_comp"] = 1

                human_id = None
                virus_id = None
                if config1["left_rev"] == 1:
                    if config1["left_id"] == "supple":
                        human_id = "read1"
                        virus_id = "read1_sa"
                    else:
                        human_id = "read1_sa"
                        virus_id = "read1"

                    res["act"][virus_id]["comp"] = 1

                    if config2["left_id"] == "supple":
                        human_id = "read2"
                        virus_id = "read2_sa"
                    else:
                        human_id = "read2_sa"
                        virus_id = "read2"

                    res["act"][virus_id]["comp"]        = 1
                    res["act"][virus_id]["strand_comp"] = 1
                    res["act"][human_id]["strand_comp"] = 1
                elif config1["right_rev"] == 1:
                    if config1["left_id"] == "supple":
                        human_id = "read1_sa"
                        virus_id = "read1"
                    else:
                        human_id = "read1"
                        virus_id = "read1_sa"

                    res["act"][virus_id]["comp"] = 1

                    if config2["left_id"] == "supple":
                        human_id = "read2_sa"
                        virus_id = "read2"
                    else:
                        human_id = "read2"
                        virus_id = "read2_sa"

                    res["act"][virus_id]["comp"]        = 1
                    res["act"][virus_id]["strand_comp"] = 1
                    res["act"][human_id]["strand_comp"] = 1
            else:
                res["left"]     = config2
                res["left_id"]  = "read2"
                res["right"]    = config1
                res["right_id"] = "read1"

                res["act"]["read1"]["strand_comp"]    = 1
                res["act"]["read1_sa"]["strand_comp"] = 1

                human_id = None
                virus_id = None
                if config1["left_rev"] == 1:
                    if config1["left_id"] == "supple":
                        human_id = "read1"
                        virus_id = "read1_sa"
                    else:
                        human_id = "read1_sa"
                        virus_id = "read1"

                    res["act"][virus_id]["comp"]        = 1
                    res["act"][virus_id]["strand_comp"] = 1
                    res["act"][human_id]["strand_comp"] = 1

                    if config2["left_id"] == "supple":
                        human_id = "read2"
                        virus_id = "read2_sa"
                    else:
                        human_id = "read2_sa"
                        virus_id = "read2"

                    res["act"][virus_id]["comp"] = 1
                elif config1["right_rev"] == 1:
                    if config1["left_id"] == "supple":
                        human_id = "read1_sa"
                        virus_id = "read1"
                    else:
                        human_id = "read1"
                        virus_id = "read1_sa"

                    res["act"][virus_id]["comp"]        = 1
                    res["act"][virus_id]["strand_comp"] = 1
                    res["act"][human_id]["strand_comp"] = 1

                    if config2["left_id"] == "supple":
                        human_id = "read2_sa"
                        virus_id = "read2"
                    else:
                        human_id = "read2"
                        virus_id = "read2_sa"

                    res["act"][virus_id]["comp"] = 1

            return res


    def toFragment(self):
        chrVirus = self.para.virus
        vlen = self.vlen()

        circular_stat = 1
        if self.para.linear == True:
            circular_stat = 0

        res = {
            # Equals config_stat;   
            "stat" : 0,
            "rnam" : "NA",
            "left" : {
                "chro"   : "NA",
                "strand" : ".",
                "start"  : 0,
                "stop"   : 0,
                "sequ"   : "NA",
                "qual"   : "NA",
                "rev"    : 0,
                "circ"   : 0,
            },
            "right" : {
                "chro"   : "NA",
                "strand" : ".",
                "start"  : 0,
                "stop"   : 0,
                "sequ"   : "NA",
                "qual"   : "NA",
                "rev"    : 0,
                "circ"   : 0,
            },
            "homolog" : {
                "length" : -1,
                "human"  : {
                    "chro"   : "NA",
                    "strand" : "NA",
                    "start"  : 0,
                    "stop"   : 0,
                },
                "virus" : {
                    "chro"   : "NA",
                    "strand" : "NA",
                    "start"  : 0,
                    "stop"   : 0,
                },
            },

            "gap" : {
                "gap_flag"   : 0,
                "gap_chro"   : "NA",
                "gap_strand" : ".",
                "gap_start"  : 0,
                "gap_stop"   : 0,
            },

            # If has breakpoint:
            # 0: blur;
            # 1: has a clear breakpoint;
            # 2: breakpoint lies in the homolog region;
            "breakpoint_stat" : 0,

            # set 1 if there are two Homolog records and
            # they are discordant;
            "homo_discord" : 0,
        }

        ###################################################
        ##             Single End Sequencing
        ###################################################

        if self.singleEndSeq():
            pair_align = self.config_OK()

            # ==================== Process the data by $act ====================
            posi = None
            read1 = "left"
            act = pair_align["act"]

            # Read1 complement
            if act["read1"]["comp"] == 1:
                if pair_align[read1]["left_id"] == "primary":
                    posi = "left"

                if pair_align[read1]["right_id"] == "primary":
                    posi = "right"

                sequ = pair_align[read1][posi]["sequ"]
                pair_align[read1][posi]["sequ"] = compSeq(sequ)

            # Read1 strand complement
            if act["read1"]["strand_comp"] == 1:
                
                #$res->{act_read1_strand_comp} = "Yes";

                if pair_align[read1]["left_id"] == "primary":
                    posi = "left"

                if pair_align[read1]["right_id"] == "primary":
                    posi = "right"

                strand = pair_align[read1][posi]["strand"]

                if strand == "+":
                    pair_align[read1][posi]["strand"] = "-"
                else:
                    pair_align[read1][posi]["strand"] = "+"

            # Read1_sa reverse
            # There's no way that ($act->{read1_sa}->{rev} == 1). The reason is same as above.

            # Read1_sa complement
            if act["read1_sa"]["comp"] == 1:
                # $res->{act_read1_strand_comp} = "Yes";

                if pair_align[read1]["left_id"] == "supple":
                    posi = "left"

                if pair_align[read1]["right_id"] == "supple":
                    posi = "right"

                sequ = pair_align[read1][posi]["sequ"]
                pair_align[read1][posi]["sequ"] = compSeq(sequ)

            # Read1_sa strand complement
            if act["read1_sa"]["strand_comp"] == 1:
                # $res->{act_read1_sa_strand_comp} = "Yes";

                if pair_align[read1]["left_id"] == "supple":
                    posi = "left"

                if pair_align[read1]["right_id"] == "supple":
                    posi = "right"

                strand = pair_align[read1][posi]["strand"]

                if strand == "+":
                    pair_align[read1][posi]["strand"] = "-"
                else:
                    pair_align[read1][posi]["strand"] = "+"

            # ====================     Processing Done!     ====================

            res["stat"] = 1
            res["breakpoint_stat"] = 1

            if pair_align["left"]["homolog"]["length"] < 0:
                pair_align["left"]["homolog"]["length"] = 0

            # Fixed in python version
            if pair_align["left"]["homolog"]["length"] > 0:
                res["breakpoint_stat"] = 2


            res["rnam"] = self.rnam()
            res["homolog"] = pair_align["left"]["homolog"]

            # Homolog was calculated by alignConfiguration in SingleAlign and its strand has yet to processed!
            if pair_align["left"]["left_spec"] == "human":
                res["homolog"]["human"]["strand"] = pair_align["left"]["left"]["strand"]
                res["homolog"]["virus"]["strand"] = pair_align["left"]["right"]["strand"]
            else:
                res["homolog"]["human"]["strand"] = pair_align["left"]["right"]["strand"]
                res["homolog"]["virus"]["strand"] = pair_align["left"]["left"]["strand"]

            # Left Fragment
            # fragL = pair_align["left"]["left"].copy()
            fragL = copy.deepcopy(pair_align["left"]["left"])
            del fragL["rnam"]
            del fragL["flag"]
            del fragL["cigar"]
            del fragL["mapq"]

            sequ_tmp  = pair_align["left"]["left"]["sequ"]
            qual_tmp  = pair_align["left"]["left"]["qual"]
            cigar_tmp = pair_align["left"]["left"]["cigar"]
            prune_tmp = pruneSequ(sequ_tmp, qual_tmp, cigar_tmp)

            fragL["sequ"] = prune_tmp["sequ"]
            fragL["qual"] = prune_tmp["qual"]
            fragL["rev"]  = pair_align["left"]["left_rev"]

            res["left"] = fragL

            # Right Fragment
            fragR = copy.deepcopy(pair_align["left"]["right"])
            del fragR["rnam"]
            del fragR["flag"]
            del fragR["cigar"]
            del fragR["mapq"]

            sequ_tmp  = pair_align["left"]["right"]["sequ"]
            qual_tmp  = pair_align["left"]["right"]["qual"]
            cigar_tmp = pair_align["left"]["right"]["cigar"]
            prune_tmp = pruneSequ(sequ_tmp, qual_tmp, cigar_tmp)

            fragR["sequ"] = prune_tmp["sequ"]
            fragR["qual"] = prune_tmp["qual"]
            fragR["rev"]  = pair_align["left"]["right_rev"]

            res["right"] = fragR

            return res

        ###################################################
        ##             Pair End Sequencing
        ###################################################

        # First process reads without supplementary alignments
        if self.is_dummy():
            res["rnam"] = self.read1.rnam

            human = None
            virus = None
            if self.read1.chro == chrVirus:
                virus = self.read1.toSingleAlignRef()
                human = self.read2.toSingleAlignRef()
            elif self.read2.chro == chrVirus:
                virus = self.read2.toSingleAlignRef()
                human = self.read1.toSingleAlignRef()

            prune_human = pruneSequ(human["sequ"], human["qual"], human["cigar"])
            prune_virus = pruneSequ(virus["sequ"], virus["qual"], virus["cigar"])

            if human["strand"] == "+":
                res["left"]["chro"]   = human["chro"]
                res["left"]["strand"] = "+"
                res["left"]["start"]  = human["start"]
                res["left"]["stop"]   = human["stop"]

                res["left"]["sequ"]  = prune_human["sequ"]
                res["left"]["qual"]  = prune_human["qual"]
                res["right"]["chro"] = chrVirus

                # Flag the Virus status:
                res["right"]["circ"] = virus["circ"]

                if virus["strand"] == "+":
                    res["right"]["strand"] = "-"
                    res["right"]["start"]  = abs(virus["stop"])
                    res["right"]["stop"]   = abs(virus["start"])

                    res["right"]["sequ"] = revCompSequ(prune_virus["sequ"])
                    res["right"]["qual"]  = prune_virus["qual"][::-1]
                    res["right"]["rev"]   = 1
                elif virus["strand"] == "-":
                    res["right"]["strand"] = "+"
                    res["right"]["start"]  = abs(virus["start"])
                    res["right"]["stop"]   = abs(virus["stop"])

                    res["right"]["sequ"] = prune_virus["sequ"]
                    res["right"]["qual"] = prune_virus["qual"]
                    res["right"]["rev"]  = 0
            elif human["strand"] == "-":
                res["right"]["chro"]   = human["chro"]
                res["right"]["strand"] = "+"
                res["right"]["start"]  = human["start"]
                res["right"]["stop"]   = human["stop"]

                res["right"]["sequ"]   = prune_human["sequ"]
                res["right"]["qual"]   = prune_human["qual"]
                res["left"]["chro"]    = chrVirus

                # Flag the Virus status:
                res["left"]["circ"] = virus["circ"]

                if virus["strand"] == "+":
                    res["left"]["strand"] = "+"
                    res["left"]["start"]  = abs(virus["start"])
                    res["left"]["stop"]   = abs(virus["stop"])

                    res["left"]["sequ"] = prune_virus["sequ"]
                    res["left"]["qual"] = prune_virus["qual"]
                    res["left"]["rev"]  = 0
                elif virus["strand"] == "-":
                    res["left"]["strand"] = "-"
                    res["left"]["start"]  = abs(virus["stop"])
                    res["left"]["stop"]   = abs(virus["start"])

                    res["left"]["sequ"] = revCompSequ(prune_virus["sequ"])
                    res["left"]["qual"] = prune_virus["qual"][::-1]
                    res["left"]["rev"]  = 1

            return res

        ############# /////////////////// ###############

        res["gap"] = {
            "gap_flag"   : 0,
            "gap_chro"   : "NA",
            "gap_strand" : ".",
            "gap_start"  : 0,
            "gap_stop"   : 0,
            "gap_rev"    : 0,
        }

        # If there's supplementary alignment(s);
        pair_align = self.config_OK();

        # Check some status first;
        if pair_align["config_stat"] == 0:
            return res

        # ==================== Process the data by $act ====================
        read1 = None
        read2 = None
        posi  = None

        if pair_align["left_id"] == "read1":
            read1 = "left"
            read2 = "right"
        else:
            read1 = "right"
            read2 = "left"

        act = pair_align["act"]

        # Read1 reverse
        if act["read1"]["rev"] == 1:
            # print("\n----- Line: 1333 -----\n")

            # $res->{act_read1_rev} = "Yes";

            if pair_align["%s_type" %(read1)] == "SingleAlignRef":
                # print("\n----- Line: 1338 -----\n")

                sta = pair_align[read1]["start"]
                sto = pair_align[read1]["stop"]
                pair_align[read1]["start"] = sto
                pair_align[read1]["stop"]  = sta

                sequ = pair_align[read1]["sequ"]
                qual = pair_align[read1]["qual"]
                pair_align[read1]["sequ"] = sequ[::-1]
                pair_align[read1]["qual"] = qual[::-1]

                cigar = pair_align[read1]["cigar"]
                cigar_rev = set_reverse_cigar(cigar)
                pair_align[read1]["cigar"] = cigar_rev

                # There's no way that ($pair_align->{"$read1\_type"} eq "config"), because all virus "reverse" operation in "config" hash map has been done by "alignConfiguration" function;

        # Read1 complement
        if act["read1"]["comp"] == 1:
            # print("\n----- Line: 1358 -----\n")

            # $res->{act_read1_comp} = "Yes";

            if pair_align["%s_type" %(read1)] == "SingleAlignRef":
                sequ = pair_align[read1]["sequ"]
                pair_align[read1]["sequ"] = compSeq(sequ)

            if pair_align["%s_type" %(read1)] == "config":
                posi = None

                if pair_align[read1]["left_id"] == "primary":
                    posi = "left"

                if pair_align[read1]["right_id"] == "primary":
                    posi = "right"

                sequ = pair_align[read1][posi]["sequ"]
                pair_align[read1][posi]["sequ"] = compSeq(sequ)

        # Read1 strand complement
        if act["read1"]["strand_comp"] == 1:
            # print("\n----- Line: 1380 -----\n")

            # $res->{act_read1_strand_comp} = "Yes";

            if pair_align["%s_type" %(read1)] == "SingleAlignRef":
                strand = pair_align[read1]["strand"]

                if strand == "+":
                    pair_align[read1]["strand"] = "-"
                else:
                    pair_align[read1]["strand"] = "+"

            if pair_align["%s_type" %(read1)] == "config":
                if pair_align[read1]["left_id"] == "primary":
                    posi = "left"

                if pair_align[read1]["right_id"] == "primary":
                    posi = "right"

                strand = pair_align[read1][posi]["strand"]

                if strand == "+":
                    pair_align[read1][posi]["strand"] = "-"
                else:
                    pair_align[read1][posi]["strand"] = "+"

        # Read1_sa reverse
        # There's no way that ($act->{read1_sa}->{rev} == 1). The reason is same as above.

        # Read1_sa complement
        if act["read1_sa"]["comp"] == 1:
            # print("\n----- Line: 1411 -----\n")

            # $res->{act_read1_sa_comp} = "Yes";

            if pair_align[read1]["left_id"] == "supple":
                posi = "left"

            if pair_align[read1]["right_id"] == "supple":
                posi = "right"

            sequ = pair_align[read1][posi]["sequ"]
            pair_align[read1][posi]["sequ"] = compSeq(sequ)

        # Read1_sa strand complement
        if act["read1_sa"]["strand_comp"] == 1:
            # print("\n----- Line: 1426 -----\n")

            # $res->{act_read1_sa_strand_comp} = "Yes";

            if pair_align[read1]["left_id"] == "supple":
                posi = "left"

            if pair_align[read1]["right_id"] == "supple":
                posi = "right"

            strand = pair_align[read1][posi]["strand"]

            if strand == "+":
                pair_align[read1][posi]["strand"] = "-"
            else:
                pair_align[read1][posi]["strand"] = "+"

        # Read2 reverse
        if act["read2"]["rev"] == 1:
            # print("\n----- Line: 1445 -----\n")

            # $res->{act_read2_rev} = "Yes";

            if pair_align["%s_type" %(read2)] == "SingleAlignRef":
                sta = pair_align[read2]["start"]
                sto = pair_align[read2]["stop"]
                pair_align[read2]["start"] = sto
                pair_align[read2]["stop"]  = sta

                sequ = pair_align[read2]["sequ"]
                qual = pair_align[read2]["qual"]
                pair_align[read2]["sequ"] = sequ[::-1]
                pair_align[read2]["qual"] = qual[::-1]

                cigar = pair_align[read2]["cigar"]
                cigar_rev = set_reverse_cigar(cigar)
                pair_align[read2]["cigar"] = cigar_rev

            # There's no way that ($pair_align->{"$read2\_type"} eq "config");

        # Read2 complement
        if act["read2"]["comp"] == 1:
            # print("\n----- Line: 1468 -----\n")

            # $res->{act_read2_comp} = "Yes";

            if pair_align["%s_type" %(read2)] == "SingleAlignRef":
                sequ = pair_align[read2]["sequ"]
                pair_align[read2]["sequ"] = compSeq(sequ)

            if pair_align["%s_type" %(read2)] == "config":
                if pair_align[read2]["left_id"] == "primary":
                    posi = "left"

                if pair_align[read2]["right_id"] == "primary":
                    posi = "right"

                sequ = pair_align[read2][posi]["sequ"]
                pair_align[read2][posi]["sequ"] = compSeq(sequ)

        # Read2 strand complement
        if act["read2"]["strand_comp"] == 1:
            # print("\n----- Line: 1488 -----\n")

            # $res->{act_read2_strand_comp} = "Yes";

            if pair_align["%s_type" %(read2)] == "SingleAlignRef":
                strand = pair_align[read2]["strand"]

                if strand == "+":
                    pair_align[read2]["strand"] = "-"
                else:
                    pair_align[read2]["strand"] = "+"

            if pair_align["%s_type" %(read2)] == "config":
                if pair_align[read2]["left_id"] == "primary":
                    posi = "left"

                if pair_align[read2]["right_id"] == "primary":
                    posi = "right"

                strand = pair_align[read2][posi]["strand"]

                if strand == "+":
                    pair_align[read2][posi]["strand"] = "-"
                else:
                    pair_align[read2][posi]["strand"] = "+"

        # Read2_sa reverse
        # There's no way that ($act->{read2_sa}->{rev} == 1);

        # Read2_sa complement
        if act["read2_sa"]["comp"] == 1:
            # print("\n----- Line: 1519 -----\n")

            # $res->{act_read2_sa_comp} = "Yes";

            if pair_align[read2]["left_id"] == "supple":
                posi = "left"

            if pair_align[read2]["right_id"] == "supple":
                posi = "right"

            sequ = pair_align[read2][posi]["sequ"]
            pair_align[read2][posi]["sequ"] = compSeq(sequ)

        # Read2_sa strand complement
        if act["read2_sa"]["strand_comp"] == 1:
            # print("\n----- Line: 1534 -----\n")

            # $res->{act_read2_sa_strand_comp} = "Yes";

            if pair_align[read2]["left_id"] == "supple":
                posi = "left"

            if pair_align[read2]["right_id"] == "supple":
                posi = "right"

            strand = pair_align[read2][posi]["strand"]

            if strand == "+":
                pair_align[read2][posi]["strand"] = "-"
            else:
                pair_align[read2][posi]["strand"] = "+"

        # ====================     Processing Done!     ====================

        res["stat"] = 1
        res["breakpoint_stat"] = 1

        #
        if pair_align["left_type"] == "config":
            # print("\n----- Line: 1558 -----\n")

            if pair_align["left"]["homolog"]["length"] < 0:
                pair_align["left"]["homolog"]["length"] = 0

            if pair_align["left"]["homolog"]["length"] > 0:
                res["breakpoint_stat"] = 2

        if pair_align["right_type"] == "config":
            # print("\n----- Line: 1564 -----\n")

            if pair_align["right"]["homolog"]["length"] < 0:
                pair_align["right"]["homolog"]["length"] = 0

            if pair_align["right"]["homolog"]["length"] > 0:
                res["breakpoint_stat"] = 2

        if pair_align["left_type"] == "config" and pair_align["right_type"] == "SingleAlignRef":
            # print("\n----- Line: 1570 -----\n")

            res["rnam"] = pair_align["right"]["rnam"]
            res["homolog"] = pair_align["left"]["homolog"]

            if res["homolog"]["length"] > 0:
                res["breakpoint_stat"] = 2

            # Homolog was calculated by alignConfiguration in SingleAlign and its strand has yet to processed!
            if(pair_align["left"]["left_spec"] == "human"):
                res["homolog"]["human"]["strand"] = pair_align["left"]["left"]["strand"]
                res["homolog"]["virus"]["strand"] = pair_align["left"]["right"]["strand"]
            else:
                res["homolog"]["human"]["strand"] = pair_align["left"]["right"]["strand"]
                res["homolog"]["virus"]["strand"] = pair_align["left"]["left"]["strand"]

            # Circular status could only exists in "SingleAlignRef";
            # $res->{right}->{circ} = $pair_align->{right}->{circ};
            # $res->{case} = 1;

            left_sequ  = pair_align["left"]["right"]["sequ"]
            left_qual  = pair_align["left"]["right"]["qual"]
            left_cigar = pair_align["left"]["right"]["cigar"]
            left_prune = pruneSequ(left_sequ, left_qual, left_cigar)

            ref_left = {
                "chro"   : pair_align["left"]["right"]["chro"],
                "strand" : pair_align["left"]["right"]["strand"],
                "start"  : pair_align["left"]["right"]["start"],
                "stop"   : pair_align["left"]["right"]["stop"],
                "sequ"   : left_prune["sequ"],
                "qual"   : left_prune["qual"],
                "rev"    : pair_align["left"]["right_rev"],

                # Left part is from a config object;    
                "circ"   : 0,   
            }

            right_sequ  = pair_align["right"]["sequ"]
            right_qual  = pair_align["right"]["qual"]
            right_cigar = pair_align["right"]["cigar"]
            right_prune = pruneSequ(right_sequ, right_qual, right_cigar)

            ref_right = {
                "chro"   : pair_align["right"]["chro"],
                "strand" : pair_align["right"]["strand"],
                "start"  : abs(pair_align["right"]["start"]),
                "stop"   : abs(pair_align["right"]["stop"]),
                "sequ"   : right_prune["sequ"],
                "qual"   : right_prune["qual"],

                # They have the same value; 
                "rev"    : pair_align["left"]["right_rev"],
                "circ"   : pair_align["right"]["circ"],
            }

            #################

            fragR, gapR = mergeSegs(ref_left, ref_right, vlen, circular_stat)
            res["right"] = fragR
            res["gap"]   = gapR

            #################

            fragL = copy.deepcopy(pair_align["left"]["left"])
            del fragL["rnam"]
            del fragL["flag"]
            del fragL["cigar"]
            del fragL["mapq"]

            sequ_tmp  = pair_align["left"]["left"]["sequ"]
            qual_tmp  = pair_align["left"]["left"]["qual"]
            cigar_tmp = pair_align["left"]["left"]["cigar"]
            prune_tmp = pruneSequ(sequ_tmp, qual_tmp, cigar_tmp)

            fragL["sequ"] = prune_tmp["sequ"]
            fragL["qual"] = prune_tmp["qual"]
            fragL["rev"]  = pair_align["left"]["left_rev"]

            # In this config, fragL's circ must be 0, whether it's from HBV or not;
            # $fragL->{circ} = 0;
            res["left"] = fragL

            return res
        elif pair_align["left_type"] == "SingleAlignRef" and pair_align["right_type"] == "config":
            # print("\n----- Line: 1652 -----\n")

            res["rnam"]    = pair_align["left"]["rnam"]
            res["homolog"] = pair_align["right"]["homolog"]

            if res["homolog"]["length"] > 0:
                res["breakpoint_stat"] = 2

            # Circular status could only exists in "SingleAlignRef";
            # Homolog was calculated by alignConfiguration in SingleAlign and its strand has yet to processed!
            if pair_align["right"]["left_spec"] == "human":
                res["homolog"]["human"]["strand"] = pair_align["right"]["left"]["strand"]
                res["homolog"]["virus"]["strand"] = pair_align["right"]["right"]["strand"]
            else:
                res["homolog"]["human"]["strand"] = pair_align["right"]["right"]["strand"]
                res["homolog"]["virus"]["strand"] = pair_align["right"]["left"]["strand"]

            left_sequ  = pair_align["left"]["sequ"]
            left_qual  = pair_align["left"]["qual"]
            left_cigar = pair_align["left"]["cigar"]
            left_prune = pruneSequ(left_sequ, left_qual, left_cigar)

            ref_left = {
                "chro"   : pair_align["left"]["chro"],
                "strand" : pair_align["left"]["strand"],
                "start"  : abs(pair_align["left"]["start"]),
                "stop"   : abs(pair_align["left"]["stop"]),
                "sequ"   : left_prune["sequ"],
                "qual"   : left_prune["qual"],
                "rev"    : pair_align["right"]["left_rev"],
                "circ"   : pair_align["left"]["circ"],
            }

            right_sequ  = pair_align["right"]["left"]["sequ"]
            right_qual  = pair_align["right"]["left"]["qual"]
            right_cigar = pair_align["right"]["left"]["cigar"]
            right_prune = pruneSequ(right_sequ, right_qual, right_cigar)

            ref_right = {
                "chro"   : pair_align["right"]["left"]["chro"],
                "strand" : pair_align["right"]["left"]["strand"],
                "start"  : pair_align["right"]["left"]["start"],
                "stop"   : pair_align["right"]["left"]["stop"],
                "sequ"   : right_prune["sequ"],
                "qual"   : right_prune["qual"],
                "rev"    : pair_align["right"]["left_rev"],
                "circ"   : 0,
            }

            ##################

            fragL, gapL = mergeSegs(ref_left, ref_right, vlen, circular_stat)
            res["left"] = fragL
            res["gap"]  = gapL

            ##################

            fragR = copy.deepcopy(pair_align["right"]["right"])
            del fragR["rnam"]
            del fragR["flag"]
            del fragR["cigar"]
            del fragR["mapq"]

            sequ_tmp  = pair_align["right"]["right"]["sequ"]
            qual_tmp  = pair_align["right"]["right"]["qual"]
            cigar_tmp = pair_align["right"]["right"]["cigar"]
            prune_tmp = pruneSequ(sequ_tmp, qual_tmp, cigar_tmp)

            fragR["sequ"] = prune_tmp["sequ"]
            fragR["qual"] = prune_tmp["qual"]
            fragR["rev"]  = pair_align["right"]["right_rev"]

            res["right"] = fragR
            return res
        elif pair_align["left_type"] == "config" and pair_align["right_type"] == "config":
            # print("\n----- Line: 1724 -----\n")

            if pair_align["left"]["homolog"]["length"] != pair_align["right"]["homolog"]["length"]:
                res["homo_discord"] = 1
                return res

            res["rnam"]          = pair_align["left"]["left"]["rnam"]
            res["left"]["chro"]  = pair_align["left"]["left"]["chro"]
            res["right"]["chro"] = pair_align["left"]["right"]["chro"]
            res["homolog"] = pair_align["left"]["homolog"]
            # res["case"] = 3

            # Homolog was calculated by alignConfiguration in SingleAlign and its strand has yet to processed!
            if pair_align["left"]["left_spec"] == "human":
                res["homolog"]["human"]["strand"] = pair_align["left"]["left"]["strand"]
                res["homolog"]["virus"]["strand"] = pair_align["left"]["right"]["strand"]
            else:
                res["homolog"]["human"]["strand"] = pair_align["left"]["right"]["strand"]
                res["homolog"]["virus"]["strand"] = pair_align["left"]["left"]["strand"]
            
            # Left fragment
            left_sequ  = pair_align["left"]["left"]["sequ"]
            left_qual  = pair_align["left"]["left"]["qual"]
            left_cigar = pair_align["left"]["left"]["cigar"]
            left_prune = pruneSequ(left_sequ, left_qual, left_cigar)

            ref_left = {
                "chro"   : pair_align["left"]["left"]["chro"],
                "strand" : pair_align["left"]["left"]["strand"],
                "start"  : pair_align["left"]["left"]["start"],
                "stop"   : pair_align["left"]["left"]["stop"],
                "sequ"   : left_prune["sequ"],
                "qual"   : left_prune["qual"],
                "rev"    : pair_align["left"]["left_rev"],
                "circ"   : 0,
            }

            right_sequ  = pair_align["right"]["left"]["sequ"]
            right_qual  = pair_align["right"]["left"]["qual"]
            right_cigar = pair_align["right"]["left"]["cigar"]
            right_prune = pruneSequ(right_sequ, right_qual, right_cigar)

            ref_right = {
                "chro"   : pair_align["right"]["left"]["chro"],
                "strand" : pair_align["right"]["left"]["strand"],
                "start"  : pair_align["right"]["left"]["start"],
                "stop"   : pair_align["right"]["left"]["stop"],
                "sequ"   : right_prune["sequ"],
                "qual"   : right_prune["qual"],
                "rev"    : pair_align["right"]["left_rev"],
                "circ"   : 0,
            }


            fragL, gapL = mergeSegs(ref_left, ref_right, vlen, circular_stat)
            res["left"] = fragL

            # Right fragment
            left_sequ  = pair_align["left"]["right"]["sequ"]
            left_qual  = pair_align["left"]["right"]["qual"]
            left_cigar = pair_align["left"]["right"]["cigar"]

            left_prune = pruneSequ(left_sequ, left_qual, left_cigar)
            ref_left = {
                "chro"   : pair_align["left"]["right"]["chro"],
                "strand" : pair_align["left"]["right"]["strand"],
                "start"  : pair_align["left"]["right"]["start"],
                "stop"   : pair_align["left"]["right"]["stop"],
                "sequ"   : left_prune["sequ"],
                "qual"   : left_prune["qual"],
                "rev"    : pair_align["left"]["right_rev"],
                "circ"   : 0,
            }

            right_sequ  = pair_align["right"]["right"]["sequ"]
            right_qual  = pair_align["right"]["right"]["qual"]
            right_cigar = pair_align["right"]["right"]["cigar"]

            right_prune = pruneSequ(right_sequ, right_qual, right_cigar)
            ref_right = {
                "chro"   : pair_align["right"]["right"]["chro"],
                "strand" : pair_align["right"]["right"]["strand"],
                "start"  : pair_align["right"]["right"]["start"],
                "stop"   : pair_align["right"]["right"]["stop"],
                "sequ"   : right_prune["sequ"],
                "qual"   : right_prune["qual"],
                "rev"    : pair_align["right"]["right_rev"],
                "circ"   : 0,
            }

            fragR, gapR = mergeSegs(ref_left, ref_right, vlen, circular_stat)
            res["right"] = fragR

            if pair_align["left"]["homolog"]["length"] == 0:
                res["breakpoint_stat"] = 1
            elif pair_align["left"]["homolog"]["length"] > 0:
                res["breakpoint_stat"] = 2

            return res


if __name__ == "__main__":
    import __ParseArgs
    import SingleSamAlign
    import SingleAlign

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

    # pp.pprint(rec)

    for read_tmp in reads:
        id = read_tmp
        print(">%s" %(id))

        pairAln = PairSamAlign(rec[read_tmp], args)
        print("\n ----- Read 1 primary -----")
        pairAln.read1.printClass()

        print("\n ----- Read 1 supplementary -----")
        pairAln.read1_sa.printClass()

        print("\n ----- Read 2 primary -----")
        pairAln.read2.printClass()

        print("\n ----- Read 2 supplementary -----")
        pairAln.read2_sa.printClass()

        print()
        print("The virus length is: %d" %(pairAln.vlen()))
        print("The read name is: %s" %(pairAln.rnam()))

        # pairAln.read1.mapq = 0
        flag = pairAln.mapqCheck()
        if flag == 0:
            print("The mapping looks good.")
        else:
            print("The mapping looks not good.")

        # pairAln.sign = 1
        if pairAln.flagCheck() == 0:
            print("The pair-end mapping looks good.")
        else:
            print("The pair-end mapping looks not good.")

        if pairAln.is_unmapped() == 0:
            print("No mapping was flagged as unmapped.")

        if pairAln.is_discordant() == 0:
            print("The pair-end mapping is concordant.")

        if pairAln.not_pair_mapped() == 0:
            print("The mapping is paired.")

        if pairAln.is_dummy() == 0:
            print("The mapping is not dummy.")

        if pairAln.singleEndSeq() == 0:
            print("The sequencing is pair-ended.")

        if pairAln.is_strand_abnormal() == 0:
            print("The strand config looks good.")

        if pairAln.alignQC() == 0:
            print("The mapping QC looks good.")
        
        config = pairAln.config_OK()
        pp.pprint(config)

        print("\n----- toFragment Test -----\n")

        frag = pairAln.toFragment()
        pp.pprint(frag)


        print("\n ---------- <END> ----------\n")

        # break

