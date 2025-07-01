#!/usr/bin/env python

import re
import copy
import math

import __Utils
import FragmentAlign


##### Auxilary functions:
#

intersect            = __Utils.intersect
mergeSegs            = __Utils.mergeSegs
distance             = __Utils.distance
distanceVirus        = __Utils.distanceVirus
distanceVirusPattern = __Utils.distanceVirusPattern
dustScore            = __Utils.dustScore
entropyScore         = __Utils.entropyScore

#
#####

class FragmentCluster:
    __slots__ = [
        "rnams", "frags", "frags_dict", "nfrags",
        "human_posi", "virus_posi",
        "orientation", "human", "virus",
        "breakpoint", "homolog", "para",
        "merge_stat",
    ]

    def __init__(self, frag = None):
        self.rnams  = []
        self.frags  = []
        self.frags_dict = {}
        self.nfrags = 0
        self.human_posi  = "NA"
        self.virus_posi  = "NA"
        self.orientation = "NA"
        self.human = {
            "chro"   : "NA",
            "strand" : ".",
            "start"  : 0,
            "stop"   : 0,
            "rev"    : 0,
            "circ"   : 0,
        }
        self.virus = {
            "chro"   : "NA",
            "strand" : ".",
            "start"  : 0,
            "stop"   : 0,
            "rev"    : 0,
            "circ"   : 0,
        }
        self.breakpoint = {
            "human" : {
                "chro"   : "NA",
                "strand" : ".",
                "cord"   : 0,
                "start"  : 0,
                "stop"   : 0,
            },
            "virus" : {
                "chro"   : "NA",
                "strand" : ".",
                "cord"   : 0,
                "start"  : 0,
                "stop"   : 0,
            },
            # -1: NA; 0: Anomalous;
            #  1: Not Anomalous; 
            "type" : -1
        }
        self.homolog = {
            "length" : -1,
            "human"  : {
                "chro"   : "NA",
                "strand" : "NA",
                "start"  : 0,
                "stop"   : 0,
            },
            "virus"  : {
                "chro"   : "NA",
                "strand" : "NA",
                "start"  : 0,
                "stop"   : 0,
            },
        }
        self.para       = ""
        self.merge_stat = 0

        if frag == None:
            return
        
        if not isinstance(frag, FragmentAlign.FragmentAlign):
            return
        
        self.rnams.append(frag.rnam)

        tmp = copy.deepcopy(frag)
        del tmp.orientation
        self.frags.append(tmp)
        self.frags_dict[frag.fragStr()] = ""
        self.nfrags = 1

        config = frag.configPosi()
        posiH, posiV = config["human"], config["virus"]
        self.human_posi  = posiH
        self.virus_posi  = posiV
        self.orientation = frag.orientation

        self.human = dict(
            chro   = frag.chro(posiH),
            strand = frag.strand(posiH),
            start  = frag.start(posiH),
            stop   = frag.stop(posiH),
            rev    = 0,
            circ   = 0,
        )
        self.virus = dict(
            chro   = frag.chro(posiV),
            strand = frag.strand(posiV),
            start  = frag.start(posiV),
            stop   = frag.stop(posiV),
            rev    = frag.revStat(posiV),
            circ   = frag.circStat(posiV),
        )

        self.breakpoint = {
            "human" : frag.figureBreakpoint("human"),
            "virus" : frag.figureBreakpoint("virus"),
        }

        if frag.breakpoint_stat == 0:
            self.breakpoint["type"] = 0
        
        if frag.breakpoint_stat == 1 or frag.breakpoint_stat == 2:
            self.breakpoint["type"] = 1
        
        self.para = frag.para
    

    def mergeStat(self, stat = None):
        if stat != None:
            self.merge_stat = stat
        
        return self.merge_stat
    

    def retrieveRnams(self, rnam = None):
        if rnam != None:
            self.rnams.append(rnam)
        
        return self.rnams 
    

    def retrieveFrags(self, frag = None):
        if frag != None:
            self.frags.append(frag)
            self.frags_dict[frag.fragStr()] = ""
        
        return self.frags
    

    def humanPart(self):
        return self.human
    

    def virusPart(self):
        return self.virus
    

    def virus_rev(self):
        return self.virusPart()["rev"]
    

    def virus_circ(self):
        return self.virusPart()["circ"]
    

    def isEmpty(self):
        if self.nfrags == 0:
            return 1
        else:
            return 0
    

    def fragOverlap(self, frag, spec = None):
        if spec == None:
            spec = "human"

        ss = "self.%s" %(spec)
        frag_tmp = eval(ss)
        
        sta_clus = frag_tmp["start"]
        sto_clus = frag_tmp["stop"]

        config    = frag.configPosi()
        posi      = config[spec]
        sta_frag  = frag.start(posi)
        sto_frag  = frag.stop(posi)
        circ_frag = frag.circStat(posi)

        if spec == "human":
            if sta_frag >= sta_clus and sta_frag <= sto_clus:
                return 1
            
            if sta_clus >= sta_frag and sta_clus <= sto_frag:
                return 1
            
            return 0
        
        if spec == "virus":
            if self.virus_rev() == 0:
                if self.para.linear == True:
                    if sta_frag >= sta_clus and sta_frag <= sto_frag:
                        return 1
                    
                    if sta_clus >= sta_frag and sta_clus <= sto_frag:
                        return 1
                else:
                    if self.virus_circ() == 0:
                       if circ_frag == 0:
                            if sta_frag >= sta_clus and sta_frag <= sto_frag:
                                return 1
                            
                            if sta_clus >= sta_frag and sta_clus <= sto_frag:
                                return 1
                       else:
                            flag = 0
                            cod = (sta_clus + sto_clus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1
                            
                            if flag == 0:
                                if sta_clus <= sto_frag:
                                    return 1
                            else:
                                if sto_clus >= sta_frag:
                                    return 1
                    else:
                        if circ_frag == 0:
                            flag = 0
                            cod = (sta_frag + sto_frag)/2

                            if self.para.vlen - cod < cod - 1:
                                flag = 1
                            
                            if flag == 0:
                                if sta_frag <= sto_clus:
                                    return 1
                            else:
                                if sto_frag >= sta_clus:
                                    return 1
                        else:
                            # All are circular:
                            return 1
            else:
                if self.para.vlen == True:
                    if sta_frag <= sta_clus and sta_frag >= sto_clus:
                        return 1
                    
                    if sta_clus <= sta_frag and sta_clus >= sto_frag:
                        return 1
                else:
                    if self.virus_circ() == 0:
                        if circ_frag == 0:
                            if sta_frag <= sta_clus and sta_frag >= sto_clus:
                                return 1
                            
                            if sta_clus <= sta_frag and sta_clus >= sto_frag:
                                return 1
                        else:
                            flag = 0
                            cod = (sta_clus + sto_clus)/2

                            if self.para.vlen - cod < cod - 1:
                                flag = 1
                            
                            if flag == 0:
                                if sto_clus <= sta_frag:
                                    return 1
                            else:
                                if sta_clus >= sto_frag:
                                    return 1
                    else:
                        if circ_frag == 0:
                            flag = 0
                            cod = (sta_frag + sto_frag)/2

                            if self.para.vlen - cod < cod - 1:
                                flag = 1
                            
                            if flag == 0:
                                if sto_frag <= sta_clus:
                                    return 1
                            else:
                                if sta_frag >= sto_clus:
                                    return 1
                        else:
                            # All are circular.
                            return 1
            
            return 0
    

    def frag_distance(self, frag, spec = None):
        if spec == None:
            spec = "human"
        
        dis = None
        if spec == "human":
            sta_clus = self.humanPart()["start"]
            sto_clus = self.humanPart()["stop"]

            ######################################################

            #   When breakpoint_stat == 0, the formula is wrong!
            
            cord_clus = (sta_clus + sto_clus)/2

            ######################################################

            sta_frag = frag.humanPart()["start"]
            sto_frag = frag.humanPart()["stop"]

            cord_frag = (sta_frag + sto_frag)/2
            dis = abs(cord_clus - cord_frag)
        elif spec == "virus":
            sta_clus = self.virusPart()["start"]
            sto_clus = self.virusPart()["stop"]

            cord_clus = None

            if self.para.linear == True:
                cord_clus = (sta_clus + sto_clus)/2
            else:
                if self.virus_circ() == 0:
                    cord_clus = (sta_clus + sto_clus)/2
                else:
                    cord_clus = (sta_clus + sto_clus + self.para.vlen)/2

                    # 2020-2-25 (Untested)
                    if cord_clus > self.para.vlen:
                        cord_clus -= self.para.vlen
            
            sta_frag = frag.virusPart()["start"]
            sto_frag = frag.virusPart()["stop"]

            cord_frag = None
            if self.para.linear == True:
                cord_frag = (sta_frag + sto_frag)/2
            else:
                if frag.circStat(frag.configPosi()["virus"]) == 0:
                    cord_frag = (sta_frag + sto_frag)/2
                else:
                    cord_frag = (sta_frag + sto_frag + self.para.vlen)/2

                    # 2020-2-25 (Untested)
                    if cord_frag > self.para.vlen:
                        cord_frag -= self.para.vlen
            
            dis = abs(cord_clus - cord_frag)

            if self.para.linear == True:
                return dis
            else:
                dis1 = self.para.vlen - cord_clus + cord_frag
                dis2 = self.para.vlen - cord_frag + cord_clus
                
                dis_tmp = dis1
                if dis1 >= dis2:
                    dis_tmp = dis2
                
                if dis_tmp < dis:
                    dis = dis_tmp
        
        return dis
    

    def mergeStatus(self, frag, spec = None):
        if spec == None:
            spec = "human"
        
        if self.fragOverlap(frag, spec) == 1:
            return 1
        
        dist = self.frag_distance(frag, spec)
        thre = distance(self.para, spec)

        if dist <= thre:
            return 1
        else:
            return 0
    

    def addToCluster(self, frag):
        if not isinstance(frag, FragmentAlign.FragmentAlign):
            return
        
        if self.para == "":
            self.para = frag.para
        
        vlen   = self.para.vlen
        config = frag.configPosi()
        posiH, posiV = config["human"], config["virus"]

        if self.isEmpty():
            self.rnams.append(frag.rnam)
            self.frags.append(frag)
            self.frags_dict[frag.fragStr()] = ""
            self.nfrags = 1

            self.human_posi = posiH
            self.virus_posi = posiV
            self.orientation = frag.orientation

            self.human = copy.deepcopy(frag.humanPart())
            del self.human["sequ"]
            del self.human["qual"]

            self.virus = copy.deepcopy(frag.virusPart())
            del self.virus["sequ"]
            del self.virus["qual"]

            self.breakpoint = {
                "human" : frag.figureBreakpoint("human"),
                "virus" : frag.figureBreakpoint("virus"),
            }

            if frag.breakpoint_stat == 0:
                self.breakpoint["type"] = 0
            
            if frag.breakpoint_stat == 1 or frag.breakpoint_stat == 2:
                self.breakpoint["type"] = 1
            
            return 0
        

        #	///// Check if the fragments are identical:
            
        # frag_vec = self.frags
        # nfrags   = self.nfrags

        dup = self.para.duplicates

        if dup == 1:
            if frag.fragStr() in self.frags_dict:
                return 0

            # for i in range(0, nfrags):
            #     frag_tmp = frag_vec[i]

            #     if not frag.sameAlignCord(frag_tmp):
            #         continue

            #     if frag.fragStr() in self.frags_dict:
            #         continue
                
            #     frag_res = frag.combineFrags(frag_tmp)
            #     self.frags[i] = frag_res

            #     return 0

        #	/////

        self.retrieveRnams(frag.rnam)
        self.retrieveFrags(frag)
        self.nfrags += 1

        #   if($spec eq "human")
        sta_frag = frag.start(self.human_posi)
        sta_clus = self.humanPart()["start"]

        if sta_frag < sta_clus:
            self.humanPart()["start"] = sta_frag
        

        sto_frag = frag.stop(self.human_posi)
        sto_clus = self.humanPart()["stop"]

        if sto_frag > sto_clus:
            self.humanPart()["stop"] = sto_frag


        #   if($spec eq "virus")
        sta_frag = frag.start(self.virus_posi)
        sto_frag = frag.stop(self.virus_posi)
        sta_clus = self.virusPart()["start"]
        sto_clus = self.virusPart()["stop"]

        if self.virus_rev() == 0:
            if self.para.linear == True:
                if sta_frag < sta_clus:
                    self.virusPart()["start"] = sta_frag

                if sto_frag > sto_clus:
                    self.virusPart()["stop"] = sto_frag
            else:
                if self.virus_circ() == 0 and frag.circStat(posiV) == 0:
                    # head: 0; tail: 1;

                    flag1 = 0
                    flag2 = 0

                    cod_frag = (sta_frag + sto_frag)/2
                    cod_clus = (sta_clus + sto_clus)/2

                    if vlen - cod_frag < cod_frag - 1:
                        flag1 = 1

                    if vlen - cod_clus < cod_clus - 1:
                        flag2 = 1

                    if flag1 == 0 and flag2 == 1:
                        
                        # Compare linear and circular distance first;
                        if abs(cod_frag - cod_clus) <= (vlen - cod_clus + 1) + cod_frag:
                            if sta_frag < sta_clus:
                                self.virusPart()["start"] = sta_frag

                            if sto_frag > sto_clus:
                                self.virusPart()["stop"] = sto_frag
                        else:
                            self.virusPart()["circ"] = 1
                            self.virusPart()["stop"] = sto_frag
                    elif flag1 == 1 and flag2 == 0:
                        
                        # Compare linear and circular distance first;
                        if abs(cod_frag - cod_clus) <= (vlen - cod_frag + 1) + cod_clus:
                            if sta_frag < sta_clus:
                                self.virusPart()["start"] = sta_frag

                            if sto_frag > sto_clus:
                                self.virusPart()["stop"] = sto_frag
                        else:
                            self.virusPart()["circ"] = 1
                            self.virusPart()["start"] = sta_frag
                    else:
                        if sta_frag < sta_clus:
                            self.virusPart()["start"] = sta_frag

                        if sto_frag > sto_clus:
                            self.virusPart()["stop"] = sto_frag
                elif self.virus_circ() == 1 and frag.circStat(posiV) == 0:
                    if vlen - sta_frag < sta_frag - 1:
                        if sta_frag < sta_clus:
                            self.virusPart()["start"] = sta_frag
                    else:
                        if sto_frag > sto_clus:
                            self.virusPart()["stop"] = sto_frag
                elif self.virus_circ() == 0 and frag.circStat(posiV) == 1:

                    self.virusPart()["circ"] = 1

                    if vlen - sta_clus < sta_clus - 1:
                        if sta_frag < sta_clus:
                            self.virusPart()["start"] = sta_frag

                        self.virusPart()["stop"] = sto_frag
                    else:
                        self.virusPart()["start"] = sta_frag

                        if sto_frag > sto_clus:
                            self.virusPart()["stop"] = sto_frag
                else:
                    if sta_frag < sta_clus:
                        self.virusPart()["start"] = sta_frag

                    if sto_frag > sto_clus:
                        self.virusPart()["stop"] = sto_frag
        elif self.virus_rev() == 1:
            if self.para.linear == True:
                if sta_frag > sta_clus:
                    self.virusPart()["start"] = sta_frag

                if sto_frag < sto_clus:
                    self.virusPart()["stop"] = sto_frag
            else:
                if self.virus_circ() == 0 and frag.circStat(posiV) == 0:
                    flag1 = 0
                    flag2 = 0

                    cod_frag = (sta_frag + sto_frag)/2
                    cod_clus = (sta_clus + sto_clus)/2

                    if vlen - cod_frag < cod_frag - 1:
                        flag1 = 1

                    if vlen - cod_clus < cod_clus - 1:
                        flag2 = 1

                    if flag1 == 0 and flag2 == 1:
                        # Compare linear and circular distance first;

                        if abs(cod_frag - cod_clus) <= (vlen - cod_clus + 1) + cod_frag:
                            if sta_frag > sta_clus:
                                self.virusPart()["start"] = sta_frag

                            if sto_frag < sto_clus:
                                self.virusPart()["stop"] = sto_frag
                        else:
                            self.virusPart()["circ"] = 1
                            self.virusPart()["start"] = sta_frag
                    elif flag1 == 1 and flag2 == 0:
                        # Compare linear and circular distance first;

                        if abs(cod_frag - cod_clus) <= (vlen - cod_frag + 1) + cod_clus:
                            if sta_frag > sta_clus:
                                self.virusPart()["start"] = sta_frag

                            if sto_frag < sto_clus:
                                self.virusPart()["stop"] = sto_frag
                        else:
                            self.virusPart()["circ"] = 1
                            self.virusPart()["stop"] = sto_frag
                    else:
                        if sta_frag > sta_clus:
                            self.virusPart()["start"] = sta_frag

                        if sto_frag < sto_clus:
                            self.virusPart()["stop"] = sto_frag
                elif self.virus_circ() == 1 and frag.circStat(posiV) == 0:
                    if vlen - sta_frag < sta_frag - 1:
                        if sto_frag < sto_clus:
                            self.virusPart()["stop"] = sto_frag
                    else:
                        if sta_frag > sta_clus:
                            self.virusPart()["start"] = sta_frag
                elif self.virus_circ() == 0 and frag.circStat(posiV) == 1:
                    self.virusPart()["circ"] = 1

                    if vlen - sta_clus < sta_clus - 1:
                        self.virusPart()["start"] = sta_frag

                        if sto_frag < sto_clus:
                            self.virusPart()["stop"] = sto_frag
                    else:
                        if sta_frag > sta_clus:
                            self.virusPart()["start"] = sta_frag

                        self.virusPart()["stop"] = sto_frag
                else:
                    if sta_frag > sta_clus:
                        self.virusPart()["start"] = sta_frag

                    if sto_frag < sto_clus:
                        self.virusPart()["stop"] = sto_frag


        if self.breakpoint["type"] == 0 and frag.breakpoint_stat == 0:
            sta_clus = self.humanPart()["start"]
            sto_clus = self.humanPart()["stop"]

            cod_clus = (sta_clus + sto_clus)/2

            self.breakpoint["human"]["cord"]  = cod_clus
            self.breakpoint["human"]["start"] = cod_clus
            self.breakpoint["human"]["stop"]  = cod_clus
        elif self.breakpoint["type"] == 0 and frag.breakpoint_stat != 0:
            sta = frag.breakpoint["human"]["start"]
            sto = frag.breakpoint["human"]["stop"]

            cord = (sta + sto)/2

            if sto == -1:
                cord = sta

            if sta == -1:
                cord = sto

            self.breakpoint["human"]["cord"]  = cord
            self.breakpoint["human"]["start"] = cord
            self.breakpoint["human"]["stop"]  = cord

            self.breakpoint["type"] = 1
        elif self.breakpoint["type"] != 0 and frag.breakpoint_stat != 0:
            cod_clus = self.breakpoint["human"]["cord"]
            cod_frag = frag.figureBreakpoint()["cord"]

            cord = (cod_clus + cod_frag)/2

            self.breakpoint["human"]["cord"]  = cord
            self.breakpoint["human"]["start"] = cord
            self.breakpoint["human"]["stop"]  = cord


    def sortFrags(self, spec = None):
        if spec == None:
            spec = "human"

        if self.isEmpty() or self.nfrags == 1:
            return self
        
        frags = sorted(
            self.frags,
            key = lambda x: x.figureBreakpoint(spec)["cord"]
        )
        self.frags = frags
        rnams = [frag_tmp.rnam for frag_tmp in frags]
        self.rnams = rnams

        return self
    

    def refineCluster(self):
        res = []

        if self.nfrags == 1:
            return [self]
        
        self.sortFrags("virus")
        frag_clus = FragmentCluster()
        nfrags = self.nfrags

        for i in range(0, nfrags):
            if frag_clus.isEmpty() or frag_clus.mergeStatus(self.frags[i], "virus") == 1:
                frag_clus.addToCluster(self.frags[i])
                continue

            if frag_clus.mergeStatus(self.frags[i], "virus") == 0:
                res.append(frag_clus)
                frag_clus = FragmentCluster()
                frag_clus.addToCluster(self.frags[i])
                continue

        res.append(frag_clus)
        return res


    def resortFrags(self, spec = None):
        if spec == None:
            spec = "human"

        if self.isEmpty() or self.nfrags == 1:
            return self

        frag_junction = []
        frag_split    = []
        frags = self.frags

        for frag_tmp in frags:
            if frag_tmp.breakpoint_stat == 0:
                frag_split.append(frag_tmp)
                continue

            if frag_tmp.breakpoint_stat != 0:
                frag_junction.append(frag_tmp)
                continue

        posi = None
        if spec == "human":
            posi = self.human_posi
        else:
            posi = self.virus_posi

        fragJ_sort = []
        fragS_sort = []
        if posi == "left":
            if len(frag_junction) >= 1:
                fragJ_sort = sorted(
                    frag_junction,
                    key = lambda x: x.figureBreakpoint(spec)["cord"],
                    reverse = True
                )

            if len(frag_split) >= 1:
                fragS_sort = sorted(
                    frag_split,
                    key = lambda x: x.figureBreakpoint(spec)["cord"],
                    reverse = True
                )
        else:
            if len(frag_junction) >= 1:
                fragJ_sort = sorted(
                    frag_junction,
                    key = lambda x: x.figureBreakpoint(spec)["cord"]
                )

            if len(frag_split) >= 1:
                fragS_sort = sorted(
                    frag_split,
                    key = lambda x: x.figureBreakpoint(spec)["cord"]
                )

        return fragJ_sort, fragS_sort


    ###################################################################

    # If $para->circular_stat == 0, then $self->virus_circ is always 0;

    ###################################################################
    def breakpointScore(self):
        bp_thre = self.para.bp_thre

        if self.isEmpty():
            return 0

        posi_human = self.human_posi
        posi_virus = self.virus_posi

        spec_left  = ""
        spec_right = ""
        if posi_human == "left":
            spec_left  = "human"
            spec_right = "virus"
        else:
            spec_left  = "virus"
            spec_right = "human"

        if self.nfrags == 1:
            breakpoint = self.frags[0].breakpoint
            sss = "Left_breakpoint: "

            # Left:
            chro_left   = breakpoint[spec_left]["chro"]
            strand_left = breakpoint[spec_left]["strand"]
            sta_left    = breakpoint[spec_left]["start"]
            sto_left    = breakpoint[spec_left]["stop"]

            cord = sta_left

            bpKey = ""
            bpStr = ""
            if self.frags[0].breakpoint_stat == 0:
                bpStr = "%s:%s:%d-NA" %(chro_left, strand_left, cord)
            else:
                bpStr = "%s:%s:%d" %(chro_left, strand_left, cord)

            sss += bpStr
            bpKey += bpStr + ";"

            # Right:
            chro_right   = breakpoint[spec_right]["chro"]
            strand_right = breakpoint[spec_right]["strand"]
            sta_right    = breakpoint[spec_right]["start"]
            sto_right    = breakpoint[spec_right]["stop"]

            cord = sto_right

            if self.frags[0].breakpoint_stat == 0:
                bpStr = "%s:%s:NA-%d" %(chro_right, strand_right, cord)
                sss  += "\tRight_breakpoint: %s" %(bpStr)
            else:
                bpStr = "%s:%s:%d" %(chro_right, strand_right, cord)
                sss  += "\tRight_breakpoint: %s" %(bpStr)

            bpKey += bpStr

            sco = None
            if self.frags[0].breakpoint_stat == 0:
                sco = "NA"
            else:
                sco = 1

            res = {}
            res[bpKey] = {
                "score" : sco,
                "bpStr" : sss,
            }

            return res

        # Cluster Size > 1:
        res = {}

        fragJ_sort, fragS_sort = self.resortFrags("human")
        fragJ_siz = len(fragJ_sort)
        fragS_siz = len(fragS_sort)

        # Left:
        chro_left   = self.breakpoint[spec_left]["chro"]
        strand_left = self.breakpoint[spec_left]["strand"]

        # Right:
        chro_right   = self.breakpoint[spec_right]["chro"]
        strand_right = self.breakpoint[spec_right]["strand"]

        #
        if fragJ_siz > 0 and fragS_siz == 0:
            for i in range(0, fragJ_siz):

                # Left:
                start_left = fragJ_sort[i].breakpoint[spec_left]["start"]
                stop_left  = fragJ_sort[i].breakpoint[spec_left]["stop"]
                bpStr_left = "%s:%s:%d-%d" %(chro_left, strand_left, start_left, stop_left)

                # Right:
                start_right = fragJ_sort[i].breakpoint[spec_right]["start"]
                stop_right  = fragJ_sort[i].breakpoint[spec_right]["stop"]
                bpStr_right = "%s:%s:%d-%d" %(chro_right, strand_right, start_right, stop_right)

                bpKey = "%s;%s" %(bpStr_left, bpStr_right)
                bpStr = "Left_breakpoint: %s\tRight_breakpoint: %s" %(bpStr_left, bpStr_right)

                if not bpKey in res:
                    res[bpKey] = dict(
                        score      = 1,
                        bpStr      = bpStr,
                        cord_left  = stop_left,
                        cord_right = start_right,
                        stat       = 0,
                    )
                else:
                    res[bpKey]["score"] += 1

            bpKeys = list(res.keys())
            n = len(bpKeys)

            if n == 1:
                return res

            bpKeys = sorted(
                bpKeys,
                key = lambda x: res[x]["score"],
                reverse = True
            )
            resu = {}

            # Need optimalization
            for j in range(0, n):
                if res[ bpKeys[j] ]["stat"] == 1:
                    continue

                bp_current = bpKeys[j]
                res[bp_current]["stat"] = 1
                resu[bp_current] = copy.deepcopy(res[bp_current])

                for i in range(0, n):
                    if res[ bpKeys[i] ]["stat"] == 1:
                        continue

                    if spec_left == "human":
                        dh = abs(res[bp_current]["cord_left"] - res[ bpKeys[i] ]["cord_left"])
                        dv = distanceVirus(
                            res[bp_current]["cord_right"],
                            res[ bpKeys[i] ]["cord_right"],
                            self.para.vlen,
                            1 - int(self.para.linear)
                        )

                        if dh + dv > bp_thre:
                            continue

                        resu[bp_current]["score"] += res[ bpKeys[i] ]["score"]
                        res[ bpKeys[i] ]["stat"] = 1
                    else:
                        dh = abs(res[bp_current]["cord_right"] - res[ bpKeys[i] ]["cord_right"])
                        dv = distanceVirus(
                            res[bp_current]["cord_left"],
                            res[ bpKeys[i] ]["cord_left"],
                            self.para.vlen,
                            1- int(self.para.linear)
                        )

                        if dh + dv > bp_thre:
                            continue

                        resu[bp_current]["score"] += res[ bpKeys[i] ]["score"]
                        res[ bpKeys[i] ]["stat"] = 1

            return resu
        elif fragJ_siz == 0 and fragS_siz > 0:
            bpoint = fragS_sort[0].breakpoint

            cord_left  = None
            cord_right = None

            # Left:
            cord_left = bpoint[spec_left]["start"]

            # Fix the initiation
            if spec_left == "virus" and self.virus_circ == 1:
                if self.virus_rev() == 0:
                    if cord_left >= self.para.vlen/2:
                        cord_left = 1

                if self.virus_rev() == 1:
                    if cord_left <= self.para.vlen/2:
                        cord_left = self.para.vlen

            # Right:
            cord_right = bpoint[spec_right]["stop"]

            # Fix the initiation
            if spec_right == "virus" and self.virus_circ() == 1:
                if self.virus_rev() == 0:
                    if cord_right <= self.para.vlen/2:
                        cord_right = self.para.vlen

                if self.virus_rev() == 1:
                    if cord_right >= self.para.vlen/2:
                        cord_right = 1

            for i in range(1, fragS_siz):
                bpoint = fragS_sort[i].breakpoint
                cord_left_tmp  = None
                cord_right_tmp = None

                # Left:
                cord_left_tmp = bpoint[spec_left]["start"]

                if spec_left == "virus":
                    if self.virus_rev() == 1:
                        if self.virus_circ() == 0:
                            if cord_left_tmp < cord_left:
                                cord_left = cord_left_tmp
                        else:
                            if cord_left_tmp < self.para.vlen/2:
                                continue

                            if cord_left_tmp < cord_left:
                                cord_left = cord_left_tmp
                    else:
                        if self.virus_circ() == 0:
                            if cord_left_tmp > cord_left:
                                cord_left = cord_left_tmp
                        else:
                            if cord_left_tmp > self.para.vlen/2:
                                continue

                            if cord_left_tmp > cord_left:
                                cord_left = cord_left_tmp
                else:
                    if cord_left_tmp > cord_left:
                        cord_left = cord_left_tmp

                # Right:
                cord_right_tmp = bpoint[spec_right]["stop"]

                if spec_right == "virus":
                    if self.virus_rev() == 1:
                        if self.virus_circ() == 0:
                            if cord_right_tmp > cord_right:
                                cord_right = cord_right_tmp
                        else:
                            if cord_right_tmp > self.para.vlen/2:
                                continue

                            if cord_right_tmp > cord_right:
                                cord_right = cord_right_tmp
                    else:
                        if self.virus_circ() == 0:
                            if cord_right_tmp < cord_right:
                                cord_right = cord_right_tmp
                        else:
                            if cord_right_tmp < self.para.vlen/2:
                                continue

                            if cord_right_tmp < cord_right:
                                cord_right = cord_right_tmp
                else:
                    if cord_right_tmp < cord_right:
                        cord_right = cord_right_tmp

            bpStr_left  = "%s:%s:%d-NA" %(chro_left, strand_left, cord_left)
            bpStr_right = "%s:%s:NA-%d" %(chro_right, strand_right, cord_right)

            sss = "Left_breakpoint: %s\tRight_breakpoint: %s" %(bpStr_left, bpStr_right)
            bpKey = "%s;%s" %(bpStr_left, bpStr_right)

            res[bpKey] = dict(
                score = "NA",
                bpStr = sss,
            )

            return res
        elif fragJ_siz > 0 and fragS_siz > 0:
            # Process junction reads first

            cod_left  = "stop"
            cod_right = "start"

            for i in range(0, fragJ_siz):
                bpoint = fragJ_sort[i].breakpoint

                # Left:
                cord_left  = bpoint[spec_left][cod_left]
                bpStr_left = "%s:%s:%d" %(chro_left, strand_left, cord_left)

                # Right:
                cord_right  = bpoint[spec_right][cod_right]
                bpStr_right = "%s:%s:%d" %(chro_right, strand_right, cord_right)

                sss = "Left_breakpoint: %s\tRight_breakpoint: %s" %(bpStr_left, bpStr_right);
                bpKey = "%s;%s" %(bpStr_left, bpStr_right)

                if not bpKey in res:
                    res[bpKey] = dict(
                        score      = 1,
                        bpStr      = sss,
                        cord_left  = cord_left,
                        cord_right = cord_right,
                        stat       = 0,
                    )
                else:
                    res[bpKey]["score"] += 1

            # Process disjoint reads then

            bpKeys = list(res.keys())
            n = len(bpKeys)

            for i in range(0, fragS_siz):
                bpoint = fragS_sort[i].breakpoint
                cord_left_tmp  = bpoint[spec_left]["start"]
                cord_right_tmp = bpoint[spec_right]["stop"]

                dvec = []

                for bpKey_tmp in bpKeys:
                    cord_left  = res[bpKey_tmp]["cord_left"]
                    cord_right = res[bpKey_tmp]["cord_right"]

                    if spec_left == "human":
                        dh = None

                        if cord_left_tmp > cord_left:
                            dh = "NA"
                        else:
                            dh = cord_left - cord_left_tmp

                        dv = None
                        pat = distanceVirusPattern(
                            cord_right,
                            cord_right_tmp,
                            self.para.vlen
                        )

                        if self.para.linear == True:
                            if self.virus_rev() == 0:
                                if cord_right_tmp < cord_right:
                                    dv = "NA"
                                else:
                                    dv = cord_right_tmp - cord_right
                            else:
                                if cord_right_tmp > cord_right:
                                    dv = "NA"
                                else:
                                    dv = abs(cord_right_tmp - cord_right)
                        else:
                            if self.virus_rev() == 0:
                                if pat == "linear":
                                    if cord_right_tmp < cord_right:
                                        dv = "NA"
                                    else:
                                        dv = cord_right_tmp - cord_right
                                else:
                                    if cord_right_tmp > cord_right:
                                        dv = "NA"
                                    else:
                                        dv = self.para.vlen - cord_right + cord_right_tmp
                            else:
                                if pat == "linear":
                                    if cord_right_tmp > cord_right:
                                        dv = "NA"
                                    else:
                                        dv = abs(cord_right_tmp - cord_right)
                                else:
                                    if cord_right_tmp < cord_right:
                                        dv = "NA"
                                    else:
                                        dv = self.para.vlen - cord_right_tmp + cord_right

                        if dh == "NA" or dv == "NA":
                            dvec.append(0)
                        else:
                            if dh + dv == 0 or dh + dv == 1:
                                dvec.append(1)
                            else:
                                dvec.append(1/math.sqrt(dh + dv))
                    else:
                        dh = None

                        if cord_right_tmp < cord_right:
                            dh = "NA"
                        else:
                            dh = cord_right_tmp - cord_right

                        dv  = None
                        pat = distanceVirusPattern(
                            cord_left,
                            cord_left_tmp,
                            self.para.vlen,
                        )

                        if self.para.linear == True:
                            if self.virus_rev() == 0:
                                if cord_left_tmp > cord_left:
                                    dv = "NA"
                                else:
                                    dv = abs(cord_left_tmp - cord_left)
                            else:
                                if cord_left_tmp < cord_left:
                                    dv = "NA"
                                else:
                                    dv = abs(cord_left_tmp - cord_left)
                        else:
                            if self.virus_rev() == 0:
                                if pat == "linear":
                                    if cord_left_tmp > cord_left:
                                        dv = "NA"
                                    else:
                                        dv = abs(cord_left_tmp - cord_left)
                                else:
                                    if cord_left_tmp < cord_left:
                                        dv = "NA"
                                    else:
                                        dv = self.para.vlen - cord_left_tmp + cord_left
                            else:
                                if pat == "linear":
                                    if cord_left_tmp < cord_left:
                                        dv = "NA"
                                    else:
                                        dv = abs(cord_left_tmp - cord_left)
                                else:
                                    if cord_left_tmp > cord_left:
                                        dv = "NA"
                                    else:
                                        dv = self.para.vlen - cord_left + cord_left_tmp

                        if dh == "NA" or dv == "NA":
                            dvec.append(0)
                        else:
                            if dh + dv == 0 or dh + dv == 1:
                                dvec.append(1)
                            else:
                                dvec.append(1/math.sqrt(dh + dv))

                ds = sum(dvec)
                if ds == 0:
                    continue
                else:
                    for i in range(0, n):
                        res[ bpKeys[i] ]["score"] += dvec[i]/ds

            if n == 1:
                return res

            bpKeys = sorted(
                bpKeys,
                key = lambda x: res[x]["score"],
                reverse = True,
            )
            resu = {}

            for j in range(0, n):
                if res[ bpKeys[j] ]["stat"] == 1:
                    continue

                bp_current = bpKeys[j]
                res[bp_current]["stat"] = 1

                resu[bp_current] = copy.deepcopy(res[bp_current])

                for i in range(0, n):
                    if res[ bpKeys[i] ]["stat"] == 1:
                        continue

                    if spec_left == "human":
                        dh = abs(res[bp_current]["cord_left"] - res[ bpKeys[i] ]["cord_left"])
                        dv = distanceVirus(
                            res[bp_current]["cord_right"],
                            res[ bpKeys[i] ]["cord_right"],
                            self.para.vlen,
                            1- int(self.para.linear),
                        )

                        if dh + dv > bp_thre:
                            continue

                        resu[bp_current]["score"] += res[ bpKeys[i] ]["score"]
                        res[ bpKeys[i] ]["stat"] = 1
                    else:
                        dh = abs(res[bp_current]["cord_right"] - res[ bpKeys[i] ]["cord_right"])
                        dv = distanceVirus(
                            res[bp_current]["cord_left"],
                            res[ bpKeys[i] ]["cord_left"],
                            self.para.vlen,
                            1 - int(self.para.linear)
                        )

                        if dh + dv > bp_thre:
                            continue

                        resu[bp_current]["score"] += res[ bpKeys[i] ]["score"]
                        res[ bpKeys[i] ]["stat"] = 1

            bpKeys = resu.keys()
            for bpKey_tmp in bpKeys:
                resu[bpKey_tmp]["score"] = round(resu[bpKey_tmp]["score"])

            return resu


    def homologCall(self):
        rec = {}
        frags = self.frags

        for frag_tmp in frags:
            # print("rnam: " + frag_tmp.rnam)

            if not frag_tmp.hasHomolog():
                continue

            homolog   = frag_tmp.homolog
            str_human = "%s:%s:%s-%s" %(
                homolog["human"]["chro"],
                homolog["human"]["strand"],
                str(homolog["human"]["start"]),
                str(homolog["human"]["stop"])
            )
            str_virus = "%s:%s:%s-%s" %(
                homolog["virus"]["chro"],
                homolog["virus"]["strand"],
                str(homolog["virus"]["start"]),
                str(homolog["virus"]["stop"])
            )
            ss = "%s; %s" %(str_human, str_virus)

            if not ss in rec:
                rec[ss] = dict(
                    homolog = homolog,
                    count   = 1,
                )

                continue
            else:
                rec[ss]["count"] += 1
                continue

        vec = list(rec.keys())

        res = dict(
            length = -1,
            human  = {
                "chro"   : "NA",
                "strand" : "NA",
                "start"  : 0,
                "stop"   : 0,
            },
            virus  = {
                "chro"   : "NA",
                "strand" : "NA",
                "start"  : 0,
                "stop"   : 0,
            },
        )

        if len(vec) == 0:
            return res

        res = rec[ vec[0] ]["homolog"]
        num = rec[ vec[0] ]["count"]

        if len(vec) > 1:
            for i in range(1, len(vec)):
                if rec[ vec[i] ]["count"] > num:
                    res = rec[ vec[i] ]["homolog"]

        self.homolog = res
        return res


    #
    # Check if two clusters could be merged:
    # 0: No; 1: Yes;
    #
    def mergeClusterStatus(self, clus):
        thre_human = self.para.bound_human
        thre_virus = self.para.bound_virus

        orient1 = self.orientation
        orient2 = clus.orientation

        if orient1 != orient2:
            return 0

        chr1 = self.humanPart()["chro"]
        chr2 = clus.humanPart()["chro"]

        if chr1 != chr2:
            return 0

        #
        # Both clusters have exact breakpoints:
        #
        if self.breakpoint["type"] == 1 and clus.breakpoint["type"] == 1:
            bp_ref1 = self.breakpointScore()
            bpKeys1 = list(bp_ref1.keys())

            bp_ref2 = clus.breakpointScore()
            bpKeys2 = list(bp_ref2.keys())

            bp_instersect = intersect(bpKeys1, bpKeys2)
            if len(bp_instersect) == 0:
                return 0
            else:
                return 1

        #
        # Neither clusters have exact breakpoints:
        #
        if self.breakpoint["type"] == 0 and clus.breakpoint["type"] == 0:
            sta1_human = self.humanPart()["start"]
            sto1_human = self.humanPart()["stop"]
            sta2_human = clus.humanPart()["start"]
            sto2_human = clus.humanPart()["stop"]

            flag_human = 0
            if sta2_human >= sta1_human and sta2_human <= sto1_human:
                flag_human = 1
            elif sta1_human >= sta2_human and sta1_human <= sto2_human:
                flag_human = 1
            else:
                if abs(sta2_human - sto1_human) <= thre_human:
                    flag_human = 1

                if abs(sta1_human - sto2_human) <= thre_human:
                    flag_human = 1

            sta1_virus = self.virusPart()["start"]
            sto1_virus = self.virusPart()["stop"]
            cod1  = (sta1_virus + sto1_virus)/2
            flag1 = 0

            if self.para.vlen - cod1 < cod1 - 1:
                flag1 = 1

            sta2_virus = clus.virusPart()["start"]
            sto2_virus = clus.virusPart()["stop"]
            cod2  = (sta2_virus + sto2_virus)/2
            flag2 = 0

            if self.para.vlen - cod2 < cod2 - 1:
                flag2 = 1

            flag_virus = 0
            if self.virus_rev() == 0:
                if self.para.linear == True:
                    if sta2_virus >= sta1_virus and sta2_virus <= sto1_virus:
                        flag_virus = 1
                    
                    if sta1_virus >= sta2_virus and sta1_virus <= sto2_virus:
                        flag_virus = 1
                    
                    if abs(sta2_virus - sto1_virus) <= thre_virus:
                        flag_virus = 1
                    
                    if abs(sta1_virus - sto2_virus) <= thre_virus:
                        flag_virus = 1
                else:
                    if self.virus_circ() == 0:
                        if clus.virus_circ() == 0:
                            if sta2_virus >= sta1_virus and sta2_virus <= sto1_virus:
                                flag_virus = 1
                            
                            if sta1_virus >= sta2_virus and sta1_virus <= sto2_virus:
                                flag_virus = 1
                            
                            if abs(sta2_virus - sto1_virus) <= thre_virus:
                                flag_virus = 1
                            
                            if abs(sta1_virus - sto2_virus) <= thre_virus:
                                flag_virus = 1
    
                            if self.para.vlen - sto2_virus + sta1_virus + 1 <= thre_virus:
                                flag_virus = 1
                        else:
                            # 0: Head; 1: Tail
                            flag = 0
                            cod  = (sta1_virus + sto1_virus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1
                            
                            if flag == 0:
                                if sta1_virus <= sto2_virus:
                                    flag_virus = 1
                                
                                if abs(sta1_virus - sto2_virus) <= thre_virus:
                                    flag_virus = 1
                            else:
                                if sto1_virus >= sta2_virus:
                                    flag_virus = 1
                                
                                if abs(sta2_virus - sto1_virus) <= thre_virus:
                                    flag_virus = 1
                    else:
                        if clus.virus_circ() == 0:
                            # 0: Head; 1: Tail
                            flag = 0
                            cod  = (sta2_virus + sto2_virus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1
                            
                            if flag == 0:
                                if sta2_virus <= sto1_virus:
                                    flag_virus = 1
                                
                                if abs(sta2_virus - sto1_virus) <= thre_virus:
                                    flag_virus = 1
                            else:
                                if sto2_virus >= sta1_virus:
                                    flag_virus = 1
                                
                                if abs(sta1_virus - sto2_virus) <= thre_virus:
                                    flag_virus = 1
                        else:
                            flag_virus = 1
            else:
                if self.para.linear == True:
                    if sta2_virus <= sta1_virus and sta2_virus >= sto1_virus:
                        flag_virus = 1

                    if sta1_virus <= sta2_virus and sta1_virus >= sto2_virus:
                        flag_virus = 1

                    if abs(sta2_virus - sto1_virus) <= thre_virus:
                        flag_virus = 1

                    if abs(sta1_virus - sto2_virus) <= thre_virus:
                        flag_virus = 1
                else:
                    if self.virus_circ() == 0:
                        if clus.virus_circ() == 0:
                            if sta2_virus <= sta1_virus and sta2_virus >= sto1_virus:
                                flag_virus = 1

                            if sta1_virus <= sta2_virus and sta1_virus >= sto2_virus:
                                flag_virus = 1

                            if abs(sta2_virus - sto1_virus) <= thre_virus:
                                flag_virus = 1

                            if abs(sta1_virus - sto2_virus) <= thre_virus:
                                flag_virus = 1

                            if self.para.vlen - sta2_virus + sto1_virus + 1 <= thre_virus:
                                flag_virus = 1

                            if self.para.vlen - sta1_virus + sto2_virus + 1 <= thre_virus:
                                flag_virus = 1
                        else:
                            # 0: Head; 1: Tail
                            flag = 0
                            cod = (sta1_virus + sto1_virus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1

                            if flag == 0:
                                if sto1_virus <= sta2_virus:
                                    flag_virus = 1

                                if abs(sto1_virus - sta2_virus) <= thre_virus:
                                    flag_virus = 1
                            else:
                                if sta1_virus >= sto2_virus:
                                    flag_virus = 1

                                if abs(sta1_virus - sto2_virus) <= thre_virus:
                                    flag_virus = 1
                    else:
                        if clus.virus_circ() == 0:
                            flag = 0
                            cod = (sta2_virus + sto2_virus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1

                            if flag == 0:
                                if sto2_virus <= sta1_virus:
                                    flag_virus = 1

                                if abs(sto2_virus - sta1_virus) <= thre_virus:
                                    flag_virus = 1
                            else:
                                if sta2_virus >= sto1_virus:
                                    flag_virus = 1

                                if abs(sto1_virus - sta2_virus) <= thre_virus:
                                    flag_virus = 1
                        else:
                            flag_virus = 1

            if flag_human == 1 and flag_virus == 1:
                return 1
            else:
                return 0


        #
        # One of them has exact breakpoint
        #
        clus0 = None
        clus1 = None

        if self.breakpoint["type"] == 1 and clus.breakpoint["type"] == 0:
            clus0 = clus
            clus1 = self
        elif self.breakpoint["type"] == 0 and clus.breakpoint["type"] == 1:
            clus0 = self
            clus1 = clus

        sta0_human = clus0.humanPart()["start"]
        sto0_human = clus0.humanPart()["stop"]
        sta1_human = clus1.humanPart()["start"]
        sto1_human = clus1.humanPart()["stop"]

        flag_human = 0
        if clus1.human_posi == "left":
            if sto0_human <= sto1_human and sta1_human - sto0_human <= thre_human:
                flag_human = 1
        else:
            if sta0_human >= sta1_human and sta0_human - sto1_human <= thre_human:
                flag_human = 1

        sta0_virus = clus0.virusPart()["start"]
        sto0_virus = clus0.virusPart()["stop"]
        sta1_virus = clus1.virusPart()["start"]
        sto1_virus = clus1.virusPart()["stop"]
        rev_stat   = clus1.virusPart()["rev"]
        circ_stat  = clus1.virusPart()["circ"]

        flag_virus = 0

        # Virus Left;
        if clus1.virus_posi == "left":
            if rev_stat == 0:
                if self.para.linear == True:
                    if sto0_virus <= sto1_virus and sta1_virus - sto0_virus <= thre_virus:
                        flag_virus = 1
                else:
                    if circ_stat == 0:
                        if clus0.virusPart()["circ"] == 0:
                            cod0 = (sta0_virus + sto0_virus)/2
                            flag0 = 0

                            if self.para.vlen - cod0 < cod0 - 1:
                                flag0 = 1

                            cod1 = (sta1_virus + sto1_virus)/2
                            flag1 = 0

                            if self.para.vlen - cod1 < cod1 - 1:
                                flag1 = 1

                            if flag0 + flag1 == 0 or flag0 + flag1 == 2:
                                if sto0_virus <= sto1_virus and sta1_virus - sto0_virus <= thre_virus:
                                    flag_virus = 1
                            elif flag0 == 1 and flag1 == 0:
                                if self.para.vlen - sto0_virus + sta1_virus <= thre_virus:
                                    flag_virus = 1
                        else:
                            cod1 = (sta1_virus + sto1_virus)/2
                            flag1 = 0
                            if self.para.vlen - cod1 < cod1 - 1:
                                flag1 = 1

                            if flag1 == 0:
                                if sto0_virus <= sto1_virus and sta1_virus - sto0_virus <= thre_virus:
                                    flag_virus = 1
                    else:
                        if clus0.virusPart()["circ"] == 0:
                            flag = 0
                            cod = (sta0_virus + sto0_virus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1

                            if flag == 0:
                                if sto0_virus <= sto1_virus:
                                    flag_virus = 1
                            else:
                                if sta1_virus - sto0_virus <= thre_virus:
                                    flag_virus = 1
                        else:
                            if sto0_virus <= sto1_virus:
                                flag_virus = 1
            else:
                if self.para.linear == True:
                    if sto0_virus >= sto1_virus and sto0_virus - sta1_virus <= thre_virus:
                        flag_virus = 1
                else:
                    if circ_stat == 0:
                        if clus0.virusPart()["circ"] == 0:
                            cod0  = (sta0_virus + sto0_virus)/2
                            flag0 = 0
                            if self.para.vlen - cod0 < cod0 - 1:
                                flag0 = 1

                            cod1 = (sta1_virus + sto1_virus)/2
                            flag1 = 0
                            if self.para.vlen - cod1 < cod1 - 1:
                                flag1 = 1

                            if flag0 + flag1 == 0 or flag0 + flag1 == 2:
                                if sto0_virus >= sto1_virus and sto0_virus - sta1_virus <= thre_virus:
                                    flag_virus = 1
                            elif flag0 == 0 and flag1 == 1:
                                if self.para.vlen - sta1_virus + sto0_virus <= thre_virus:
                                    flag_virus = 1
                        else:
                            cod1 = (sta1_virus + sto1_virus)/2
                            flag1 = 0
                            if self.para.vlen - cod1 < cod1 - 1:
                                flag1 = 1

                            if flag1 == 1:
                                if sto0_virus >= sto1_virus and sto0_virus - sta1_virus <= thre_virus:
                                    flag_virus = 1
                    else:
                        if clus0.virusPart()["circ"] == 0:
                            flag = 0
                            cod  = (sta0_virus + sto0_virus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1

                            if flag == 1:
                                if sto0_virus >= sto1_virus:
                                    flag_virus = 1
                            else:
                                if sto0_virus - sta1_virus <= thre_virus:
                                    flag_virus = 1
                        else:
                            if sto0_virus >= sto1_virus:
                                flag_virus = 1

        # Virus Right:
        else:
            if rev_stat == 0:
                if self.para.linear == True:
                    if sta0_virus >= sta1_virus and sta0_virus - sto1_virus <= thre_virus:
                        flag_virus = 1
                else:
                    if circ_stat == 0:
                        if clus0.virusPart()["circ"] == 0:
                            cod0 = (sta0_virus + sto0_virus)/2
                            flag0 = 0
                            if self.para.vlen -cod0 < cod0 - 1:
                                flag0 = 1

                            cod1 = (sta1_virus + sto1_virus)/2
                            flag1 = 0
                            if self.para.vlen - cod1 < cod1 - 1:
                                flag1 = 1

                            if flag0 + flag1 == 0 or flag0 + flag1 == 2:
                                if sta0_virus >= sta1_virus and sta0_virus - sto1_virus <= thre_virus:
                                    flag_virus = 1
                            elif flag0 == 0 and flag1 == 1:
                                if self.para.vlen - sto1_virus + sta0_virus <= thre_virus:
                                    flag_virus = 1
                        else:
                            cod1  = (sta1_virus + sto1_virus)/2
                            flag1 = 0
                            if self.para.vlen - cod1 < cod1 - 1:
                                flag1 = 1

                            if flag1 == 1:
                                if sta0_virus >= sta1_virus and sta0_virus - sto1_virus <= thre_virus:
                                    flag_virus = 1
                    else:
                        if clus0.virusPart()["circ"] == 0:
                            flag = 0
                            cod  = (sta0_virus + sto0_virus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1

                            if flag == 1:
                                if sta0_virus >= sta1_virus:
                                    flag_virus = 1
                            else:
                                if sta0_virus - sto1_virus <= thre_virus:
                                    flag_virus = 1
                        else:
                            if sta0_virus >= sta1_virus:
                                flag_virus = 1
            else:
                if self.para.linear == True:
                    if sta0_virus <= sta1_virus and sto1_virus - sta0_virus <= thre_virus:
                        flag_virus = 1
                else:
                    if circ_stat == 0:
                        if clus0.virusPart()["circ"] == 0:
                            cod0 = (sta0_virus + sto0_virus)/2
                            flag0 = 0
                            if self.para.vlen - cod0 < cod0 - 1:
                                flag0 = 1

                            cod1 = (sta1_virus + sto1_virus)/2
                            flag1 = 0
                            if self.para.vlen - cod1 < cod1 - 1:
                                flag1 = 1

                            if flag0 + flag1 == 0 or flag0 + flag1 == 2:
                                if sta0_virus <= sta1_virus and sto1_virus - sta0_virus:
                                    flag_virus = 1
                            else:
                                if flag0 == 1 and flag1 == 0:
                                    if self.para.vlen - sta0_virus + sto1_virus <= thre_virus:
                                        flag_virus = 1
                        else:
                            cod1 = (sta1_virus + sto1_virus)/2
                            flag1 = 0
                            if self.para.vlen - cod1 < cod1 - 1:
                                flag1 = 1

                            if flag1 == 0:
                                if sta0_virus <= sta1_virus and sto1_virus - sta0_virus <= thre_virus:
                                    flag_virus = 1
                    else:
                        if clus0.virusPart()["circ"] == 0:
                            flag = 0
                            cod  = (sta0_virus + sto0_virus)/2
                            if self.para.vlen - cod < cod - 1:
                                flag = 1

                            if flag == 0:
                                if sta0_virus <= sta1_virus:
                                    flag_virus = 1
                            else:
                                if sta0_virus - sto1_virus <= thre_virus:
                                    flag_virus = 1
                        else:
                            if sta0_virus <= sta1_virus:
                                flag_virus = 1
        
        if flag_human == 1 and flag_virus == 1:
            return 1
        else:
            return 0


    def mergeCluster(self, clus):
        import pprint
        pp = pprint.PrettyPrinter(indent = 4)


        if clus.isEmpty():
            return self

        res = FragmentCluster()

        res.para  = self.para
        rnams     = copy.deepcopy(self.rnams)
        rnams.extend(clus.rnams)
        res.rnams = rnams
        res.nfrags = self.nfrags + clus.nfrags
        frags     = copy.deepcopy(self.frags)
        frags.extend(clus.frags)
        res.frags = frags
        res.frags_dict.update(clus.frags_dict)
        res.orientation = self.orientation
        res.human_posi  = self.human_posi
        res.virus_posi  = self.virus_posi

        res.humanPart()["chro"]   = self.humanPart()["chro"]
        res.humanPart()["strand"] = self.humanPart()["strand"]
        res.breakpoint["human"]["chro"]   = self.humanPart()["chro"]
        res.breakpoint["human"]["strand"] = self.humanPart()["strand"]

        res.virusPart()["chro"]   = self.virusPart()["chro"]
        res.virusPart()["strand"] = self.virusPart()["strand"]
        res.breakpoint["virus"]["chro"]   = self.virusPart()["chro"]
        res.breakpoint["virus"]["strand"] = self.virusPart()["strand"]

        res.breakpoint["type"] = 0
        if self.breakpoint["type"] == 1 or clus.breakpoint["type"] == 1:
            res.breakpoint["type"] = 1

        # Human Coordinates:
        sta_self = self.humanPart()["start"]
        sto_self = self.humanPart()["stop"]
        sta_clus = clus.humanPart()["start"]
        sto_clus = clus.humanPart()["stop"]

        sta = sta_self
        if sta_clus < sta_self:
            sta = sta_clus

        res.humanPart()["start"] = sta

        sto = sto_self
        if sto_clus > sto_self:
            sto = sto_clus

        res.humanPart()["stop"] = sto

        # Virus Coordinates:
        res.virusPart()["rev"] = self.virus_rev()

        sta_self = self.virusPart()["start"]
        sto_self = self.virusPart()["stop"]
        sta_clus = clus.virusPart()["start"]
        sto_clus = clus.virusPart()["stop"]

        cod_self = (sta_self + sto_self)/2
        cod_clus = (sta_clus + sto_clus)/2

        # 0: Head; 1: Tail;
        flag_self = 0
        if res.para.vlen - cod_self < cod_self - 1:
            flag_self = 1

        flag_clus = 0
        if res.para.vlen - cod_clus < cod_clus - 1:
            flag_clus = 1

        if self.para.linear == False:
            if self.virus_circ() == 1 or clus.virus_circ() == 1:
                res.virusPart()["circ"] = 1

        if res.virus_rev() == 0:
            if self.para.linear == True:
                sta = sta_self
                if sta_clus < sta_self:
                    sta = sta_clus

                res.virusPart()["start"] = sta

                sto = sto_self
                if sto_clus > sto_self:
                    sto = sto_clus

                res.virusPart()["stop"] = sto
            else:
                if self.virus_circ() == 0 and clus.virus_circ() == 0:
                    if flag_self + flag_clus == 0 or flag_self + flag_clus == 2:
                        sta = sta_self
                        if sta_clus < sta_self:
                            sta = sta_clus

                        res.virusPart()["start"] = sta

                        sto = sto_self
                        if sto_clus > sto_self:
                            sto = sto_clus

                        res.virusPart()["stop"] = sto
                    else:
                        if flag_self == 1 and flag_clus == 0:
                            if res.para.vlen - cod_self + cod_clus < cod_self - cod_clus:
                                res.virusPart()["circ"] = 1

                                res.virusPart()["start"] = sta_self
                                res.virusPart()["stop"]  = sto_clus
                            else:
                                sta = sta_clus
                                if sta_self < sta_clus:
                                    sta = sta_self

                                res.virusPart()["start"] = sta

                                sto = sto_clus
                                if sto_self > sto_clus:
                                    sto = sto_self

                                res.virusPart()["stop"] = sto
                        elif flag_self == 0 and flag_clus == 1:
                            if res.para.vlen - cod_clus + cod_self < cod_clus - cod_self:
                                res.virusPart()["circ"] = 1

                                res.virusPart()["start"] = sta_clus
                                res.virusPart()["stop"]  = sto_self
                            else:
                                sta = sta_clus
                                if sta_self < sta_clus:
                                    sta = sta_self

                                res.virusPart()["start"] = sta

                                sto = sto_clus
                                if sto_self > sto_clus:
                                    sto = sto_self

                                res.virusPart()["stop"] = sto
                elif self.virus_circ() == 1 and clus.virus_circ() == 0:
                    res.virusPart()["circ"] = 1

                    if flag_clus == 0:
                        res.virusPart()["start"] = sta_self
                        res.virusPart()["stop"]  = sto_self

                        if sto_clus > sto_self:
                            res.virusPart()["stop"] = sto_clus
                    else:
                        res.virusPart()["start"] = sta_self

                        if sta_clus < sta_self:
                            res.virusPart()["start"] = sta_clus

                        res.virusPart()["stop"] = sto_self
                elif self.virus_circ() == 0 and clus.virus_circ() == 1:
                    res.virusPart()["circ"] = 1

                    if flag_self == 0:
                        res.virusPart()["start"] = sta_clus
                        res.virusPart()["stop"]  = sto_self

                        if sto_clus > sto_self:
                            res.virusPart()["stop"] = sto_clus
                    else:
                        res.virusPart()["start"] = sta_self

                        if sta_clus < sta_self:
                            res.virusPart()["start"] = sta_clus

                        res.virusPart()["stop"] = sto_clus
                else:
                    res.virusPart()["circ"] = 1

                    res.virusPart()["start"] = sta_self
                    if sta_clus < sta_self:
                        res.virusPart()["start"] = sta_clus

                    res.virusPart()["stop"] = sto_self
                    if sta_clus > sto_clus:
                        res.virusPart()["stop"] = sto_clus
        else:
            if self.para.linear == True:
                sta = sta_self
                if sta_clus > sta_self:
                    sta = sta_clus

                res.virusPart()["start"] = sta

                sto = sto_self
                if sto_clus < sto_self:
                    sto = sto_clus

                res.virusPart()["stop"] = sto
            else:
                if self.virus_circ() == 0 and clus.virus_circ() == 0:
                    if flag_self + flag_clus == 0 or flag_self + flag_clus == 2:
                        sta = sta_self
                        if sta_clus > sta_self:
                            sta = sta_clus

                        res.virusPart()["start"] = sta

                        sto = sto_self
                        if sto_clus < sto_self:
                            sto = sto_clus

                        res.virusPart()["stop"] = sto
                    else:
                        if flag_self == 1 and flag_clus == 0:
                            if res.para.vlen - cod_self + cod_clus < cod_self - cod_clus:
                                res.virusPart()["circ"] = 1

                                res.virusPart()["start"] = sta_clus
                                res.virusPart()["stop"]  = sto_self
                            else:
                                sta = sta_self
                                if sta_clus > sta_self:
                                    sta = sta_clus

                                res.virusPart()["start"] = sta

                                sto = sto_self
                                if sto_clus < sto_self:
                                    sto = sto_clus

                                res.virusPart()["stop"] = sto
                        elif flag_self == 0 and flag_clus == 1:
                            if res.para.vlen - cod_clus + cod_self < cod_clus - cod_self:
                                res.virusPart()["circ"] = 1

                                res.virusPart()["start"] = sta_self
                                res.virusPart()["stop"]  = sto_clus
                            else:
                                sta = sta_self
                                if sta_clus > sta_self:
                                    sta = sta_clus

                                res.virusPart()["start"] = sta

                                sto = sto_self
                                if sto_clus < sto_self:
                                    sto = sto_clus

                                res.virusPart()["stop"] = sto
                elif self.virus_circ() == 1 and clus.virus_circ() == 0:
                    res.virusPart()["circ"] = 1

                    if flag_clus == 0:
                        res.virusPart()["start"] = sta_self
                        if sta_clus > sta_self:
                            res.virusPart()["start"] = sta_clus

                        res.virusPart()["stop"] = sto_self
                    else:
                        res.virusPart()["start"] = sta_self
                        res.virusPart()["stop"]  = sto_self
                        if sto_clus > sto_self:
                            res.virusPart()["stop"] = sto_clus
                elif self.virus_circ() == 0 and clus.virus_circ() == 1:
                    res.virusPart()["circ"] = 1

                    if flag_self == 0:
                        res.virusPart()["start"] = sta_self
                        if sta_clus > sta_self:
                            res.virusPart()["start"] = sta_clus

                        res.virusPart()["stop"] = sto_clus
                    else:
                        res.virusPart()["start"] = sta_clus
                        res.virusPart()["stop"] = sto_self
                        if sto_clus < sto_self:
                            res.virusPart()["stop"] = sto_clus
                else:
                    res.virusPart()["circ"] = 1

                    res.virusPart()["start"] = sta_self
                    if sta_clus > sta_self:
                        res.virusPart()["start"] = sta_clus

                    res.virusPart()["stop"] = sto_self
                    if sto_clus < sto_self:
                        res.virusPart()["stop"] = sto_clus

        return res


    def mergeClusterFrags(self):
        res = dict(
            left  = "",
            right = "",
        )

        circular_stat = 1 - int(self.para.linear)

        if self.nfrags == 1:
            frag = self.frags[0]
            res["left"]  = frag.left
            res["right"] = frag.right

            return res

        frags = self.frags
        frag_merge = frags[0]
        for i in range(1, len(frags)):
            frag_merge = frag_merge.mergeFrags(frags[i], circular_stat)

        res["left"]  = frag_merge.left
        res["right"] = frag_merge.right
        
        return res


    def printFragCluster(self, filename):
        if self.isEmpty():
            return 0

        f = open(filename, "a")

        left_clus  = None
        right_clus = None
        if self.human_posi == "left":
            left_clus  = "human"
            right_clus = "virus"
        else:
            left_clus  = "virus"
            right_clus = "human"

        vlen = self.para.vlen

        ss = f"self.{left_clus}"
        clus_left = eval(ss)

        chro_left   = clus_left["chro"]
        strand_left = clus_left["strand"]
        sta_left    = clus_left["start"]
        sto_left    = clus_left["stop"]

        ss = f"self.{right_clus}"
        clus_right = eval(ss)

        chro_right   = clus_right["chro"]
        strand_right = clus_right["strand"]
        sta_right    = clus_right["start"]
        sto_right    = clus_right["stop"]

        nfrags = self.nfrags
        type   = self.breakpoint["type"]

        ss = f">Left: {chro_left}:{strand_left}:{sta_left}-{sto_left}\tRight: {chro_right}:{strand_right}:{sta_right}-{sto_right}\tCluster_Size: {nfrags}"
        f.write(ss + "\n")

        bp_ref = self.breakpointScore()
        bpKeys = list(bp_ref.keys())

        if len(bpKeys) == 1:
            bpStr = bp_ref[ bpKeys[0] ]["bpStr"]
            sco   = bp_ref[ bpKeys[0] ]["score"]

            ss = f"-{bpStr}\tScore: {sco}\tBreakpoint_type: {type}"
            f.write(ss + "\n")
        else:
            n = len(bpKeys)
            bpKeys = sorted(
                bpKeys,
                key = lambda x: bp_ref[x]["score"],
                reverse = True
            )

            for bpKey_tmp in bpKeys:
                bpStr = bp_ref[bpKey_tmp]["bpStr"]
                sco   = bp_ref[bpKey_tmp]["score"]

                ss = f"-{bpStr}\tScore: {sco}\tBreakpoint_type: {type}"
                f.write(ss + "\n")
        
        homolog = self.homologCall()
        start_left = homolog[left_clus]["start"]
        stop_left  = homolog[left_clus]["stop"]

        start_right = homolog[right_clus]["start"]
        stop_right  = homolog[right_clus]["stop"]

        if homolog["length"] > 0:
            ss = f"-Microhomolog_region: {chro_left}:{strand_left}:{start_left}-{stop_left}\t{chro_right}:{strand_right}:{start_right}-{stop_right}"
            f.write(ss + "\n")
        else:
            f.write("-Microhomolog_region: NA\n")
        
        # Merge fragments first;
        frag_sequ = self.mergeClusterFrags()
        #

        if self.para.method == 0:
            mcString = "Dust"
        else:
            mcString = "Entropy"
        
        seq_left = frag_sequ["left"]["sequ"]
        seqc = f"Left {mcString}:"

        pat = re.compile("-")
        if re.search(pat, seq_left):
            tmp = re.split(re.compile("\-+"), seq_left)
            nseg = len(tmp)

            for i in range(0, nseg):
                j = i + 1
                ss = f"Segment_{j}: "

                sco = None
                if mcString == "Dust":
                    sco = dustScore(tmp[i])
                
                if mcString == "Entropy":
                    sco = entropyScore(tmp[i])
                
                ss += str(sco)
                seqc += f" {ss};"
        else:
            sco = None
            if mcString == "Dust":
                sco = dustScore(seq_left)
            
            if mcString == "Entropy":
                sco = entropyScore(seq_left)
            
            seqc += f" Segment_1: {sco}"
        
        f.write(f"-{seqc}\n")

        seq_right = frag_sequ["right"]["sequ"]
        seqc      = f"Right {mcString}:"

        if re.search(re.compile("-"), seq_right) != None:
            tmp = re.split(re.compile("\-+"), seq_right)
            nseg = len(tmp)

            for i in range(0, nseg):
                j  = i + 1
                ss = f"Segment_{j}: "

                sco = None
                if mcString == "Dust":
                    sco = dustScore(tmp[i])
                
                if mcString == "Entropy":
                    sco = entropyScore(tmp[i])
                
                ss   += str(sco)
                seqc += f" {ss};"
        else:
            sco = None
            if mcString == "Dust":
                sco = dustScore(seq_right)
            
            if mcString == "Entropy":
                sco = entropyScore(seq_right)
            
            seqc += f" Segment_1: {sco}"
        
        f.write(f"-{seqc}\n")

        f.write("<Information Section>\n")
        for idx in range(0, self.nfrags):
            frag = self.frags[idx]
            ss   = frag.strFragmentAlign()
            f.write(f"{ss}\n")
        
        f.write("\n<Sequence Section>\n")

        # rnam width:
        lens = map(len, self.rnams)
        rwid = max(lens)

        len_ref = len("Reference_sequence")
        if len_ref > rwid:
            rwid = len_ref
        
        rwid_ss = " " * rwid

        # Sequence Width;
        swid_left  = abs(sto_left - sta_left) + 1
        swid_right = abs(sto_right - sta_right) + 1

        if clus_left["circ"] == 1:
            if self.virus_rev() == 0:
                swid_left = (vlen - sta_left + 1) + sto_left
            else:
                swid_left = (vlen - sto_left + 1) + sta_left
        
        if clus_right["circ"] == 1:
            if self.virus_rev() == 0:
                swid_right = (vlen - sta_right + 1) + sto_right
            else:
                swid_right = (vlen - sto_right + 1) + sta_right
        
        # Specie Infomation:
        left_ss  = " " * swid_left
        left_mid = round(swid_left/2 - 1)
        # ss = left_clus.upper()

        ss_left  = left_ss[0:left_mid]
        n_left   = len(left_clus)
        ss_right = left_ss[left_mid+n_left:]
        left_ss  = ss_left + left_clus.upper() + ss_right

        # ss_rep = left_ss[left_mid:left_mid+5]
        # left_ss = left_ss.replace(ss_rep, ss, 1)

        right_ss  = " " * swid_right
        right_mid = round(swid_right/2 - 1)
        ss = right_clus.upper()

        ss_left  = right_ss[0:right_mid]
        n_right  = len(right_clus)
        ss_right = right_ss[right_mid+n_right:]
        right_ss = ss_left + right_clus.upper() + ss_right

        # ss_rep = right_ss[right_mid:right_mid+5]
        # right_ss = right_ss.replace(ss_rep, ss, 1)

        f.write(f"{rwid_ss}\t{left_ss}\t{right_ss}\n")

        # Coordinate Information:
        left_ss = " "*swid_left
        ss_mid = left_ss[len(str(sta_left)):-len(str(sto_left))]
        left_ss = str(sta_left) + ss_mid + str(sto_left)

        # ss = left_ss[0:len(str(sta_left))]
        # left_ss = left_ss.replace(ss, str(sta_left), 1)

        right_ss = " "*swid_right
        ss_mid   = right_ss[len(str(sta_right)):-len(str(sto_right))]
        right_ss = str(sta_right) + ss_mid + str(sto_right)

        # idx_sta = -len(str(sto_left))
        # idx_sto = idx_sta + len(str(sto_right))
        # ss = left_ss[idx_sta:idx_sto]
        # right_ss = right_ss.replace(ss, str(sto_right), 1)

        f.write(f"{rwid_ss}\t{left_ss}\t{right_ss}\n")

        left_ss  = "|"
        left_ss += "-"*(swid_left - 2)
        left_ss += "|"

        right_ss  = "|"
        right_ss += "-"*(swid_right - 2)
        right_ss += "|"

        f.write(f"{rwid_ss}\t{left_ss}\t{right_ss}\n")

        # Print reference sequence:
        nam_header = "Reference_sequence"
        ss_header  = f'{nam_header: <{rwid}}'
        sss = "%s\t%s\t%s" %(ss_header, frag_sequ["left"]["sequ"], frag_sequ["right"]["sequ"])
        f.write(f"{sss}\n")

        for i in range(0, self.nfrags):
            frag = self.frags[i]

            # Left
            left_ss   = " " * swid_left
            flag_left = 0

            if left_clus == "human":
                idx_sta = abs(frag.start("left") - sta_left)
                idx_sto = idx_sta + len(frag.sequ("left"))
                # ss = left_ss[idx_sta:idx_sto]
                # left_ss = left_ss.replace(ss, frag.sequ("left"), 1)

                ss_left  = left_ss[0:idx_sta]
                ss_right = left_ss[idx_sto:]
                left_ss  = ss_left + frag.sequ("left") + ss_right
            # Virus:
            else:
                # cluster circ must be 1;
                if frag.circStat("left") == 1:
                    idx_sta = abs(frag.start("left") - sta_left)
                    idx_sto = idx_sta + len(frag.sequ("left"))
                    # ss = left_ss[idx_sta:idx_sto]
                    # left_ss = left_ss.replace(ss, frag.sequ("left"), 1)

                    ss_left  = left_ss[0:idx_sta]
                    ss_right = left_ss[idx_sto:]
                    left_ss  = ss_left + frag.sequ("left") + ss_right
                else:
                    if clus_left["circ"] == 0:
                        idx_sta = abs(frag.start("left") - sta_left)
                        idx_sto = idx_sta + len(frag.sequ("left"))
                        # ss = left_ss[idx_sta:idx_sto]
                        # left_ss = left_ss.replace(ss, frag.sequ("left"), 1)

                        ss_left  = left_ss[0:idx_sta]
                        ss_right = left_ss[idx_sto:]
                        left_ss  = ss_left + frag.sequ("left") + ss_right
                    else:
                        cod_tmp = (frag.start("left") + frag.stop("left"))/2
                        sta_clus = sta_left
                        sto_clus = sto_left

                        if self.virus_rev() == 0:
                            if sta_clus <= sto_clus:
                                sto_clus = sta_clus - 1

                            if frag.start("left") >= sta_clus:
                                flag_left = 1
                        else:
                            if sta_clus >= sto_clus:
                                sta_clus = sto_clus - 1

                            if frag.stop("left") >= sto_clus:
                                flag_left = 1

                        if flag_left == 0:
                            # head
                            if frag.revStat("left") == 1:
                                idx_sta = abs(frag.start("left") - sta_left)
                                idx_sto = idx_sta + len(frag.sequ("left"))
                                # ss = left_ss[idx_sta:idx_sto]
                                # left_ss.replace(ss, frag.sequ("left"), 1)

                                ss_left  = left_ss[0:idx_sta]
                                ss_right = left_ss[idx_sto:]
                                left_ss  = ss_left + frag.sequ("left") + ss_right
                            else:
                                idx_sta = vlen - sta_left + frag.start("left")
                                idx_sto = idx_sta + len(frag.sequ("left"))
                                # ss = left_ss[idx_sta:idx_sto]
                                # left_ss = left_ss.replace(ss, frag.sequ("left"), 1)

                                ss_left  = left_ss[0:idx_sta]
                                ss_right = left_ss[idx_sto:]
                                left_ss  = ss_left + frag.sequ("left") + ss_right
                        else:
                            # tail
                            if frag.revStat("left") == 0:
                                idx_sta = frag.start("left") - sta_left
                                idx_sto = idx_sta + len(frag.sequ("left"))
                                # ss = left_ss[idx_sta:idx_sto]
                                # left_ss = left_ss.replace(ss, frag.sequ("left"), 1)

                                ss_left  = left_ss[0:idx_sta]
                                ss_right = left_ss[idx_sto:]
                                left_ss  = ss_left + frag.sequ("left") + ss_right
                            else:
                                idx_sta = vlen - frag.start("left") + sta_left
                                idx_sto = idx_sta + len(frag.sequ("left"))
                                # ss = left_ss[idx_sta:idx_sto]
                                # left_ss = left_ss.replace(ss, frag.sequ("left"), 1)

                                ss_left  = left_ss[0:idx_sta]
                                ss_right = left_ss[idx_sto:]
                                left_ss  = ss_left + frag.sequ("left") + ss_right

            # Right
            right_ss = " " * swid_right
            flag_right = 0

            if right_clus == "human":
                idx_sta = abs(frag.start("right") - sta_right)
                idx_sto = idx_sta + len(frag.sequ("right"))
                # ss = right_ss[idx_sta:idx_sto]
                # right_ss = right_ss.replace(ss, frag.sequ("right"), 1)

                ss_left  = right_ss[0:idx_sta]
                ss_right = right_ss[idx_sto:]
                right_ss = ss_left + frag.sequ("right") + ss_right
            else:
                if frag.circStat("right") == 1:
                    idx_sta = abs(frag.start("right") - sta_right)
                    idx_sto = idx_sta + len(frag.sequ("right"))
                    # ss = right_ss[idx_sta:idx_sto]
                    # right_ss = right_ss.replace(ss, frag.sequ("right"), 1)

                    ss_left  = right_ss[0:idx_sta]
                    ss_right = right_ss[idx_sto:]
                    right_ss = ss_left + frag.sequ("right") + ss_right
                else:
                    if clus_right["circ"] == 0:
                        idx_sta = abs(frag.start("right") - sta_right)
                        idx_sto = idx_sta + len(frag.sequ("right"))
                        # ss = right_ss[idx_sta:idx_sto]
                        # right_ss = right_ss.replace(ss, frag.sequ("right"), 1)

                        ss_left  = right_ss[0:idx_sta]
                        ss_right = right_ss[idx_sto:]
                        right_ss = ss_left + frag.sequ("right") + ss_right
                    else:
                        cod_tmp = (frag.start("right") + frag.stop("right"))/2
                        sta_clus = sta_right
                        sto_clus = sto_right

                        if self.virus_rev() == 0:
                            if sto_clus >= sta_clus:
                                sto_clus = sta_clus - 1

                            if frag.start("right") >= sta_clus:
                                flag_right = 1
                        else:
                            if sta_clus >= sto_clus:
                                sta_clus = sto_clus - 1

                            if frag.stop("right") >= sto_clus:
                                flag_right = 1

                        if flag_right == 0:
                            # head
                            if frag.revStat("right") == 1:
                                idx_sta = abs(frag.start("right") - sta_right)
                                idx_sto = idx_sta + len(frag.sequ("right"))
                                # ss = right_ss[idx_sta:idx_sto]
                                # right_ss = right_ss.replace(ss, frag.sequ("right"), 1)

                                ss_left  = right_ss[0:idx_sta]
                                ss_right = right_ss[idx_sto:]
                                right_ss = ss_left + frag.sequ("right") + ss_right
                            else:
                                idx_sta = vlen - sta_right + frag.start("right")
                                idx_sto = idx_sta + len(frag.sequ("right"))
                                # ss = right_ss[idx_sta:idx_sto]
                                # right_ss = right_ss.replace(ss, frag.sequ("right"), 1)

                                ss_left  = right_ss[0:idx_sta]
                                ss_right = right_ss[idx_sto:]
                                right_ss = ss_left + frag.sequ("right") + ss_right
                        else:
                            # tail
                            if frag.revStat("right") == 0:
                                idx_sta = frag.start("right") - sta_right
                                idx_sto = idx_sta + len(frag.sequ("right"))
                                # ss = right_ss[idx_sta:idx_sto]
                                # right_ss = right_ss.replace(ss, frag.sequ("right"), 1)

                                ss_left  = right_ss[0:idx_sta]
                                ss_right = right_ss[idx_sto:]
                                right_ss = ss_left + frag.sequ("right") + ss_right
                            else:
                                idx_sta = vlen - frag.start("right") + sta_right
                                idx_sto = idx_sta + len(frag.sequ("right"))
                                # ss = right_ss[idx_sta:idx_sto]
                                # right_ss = right_ss.replace(ss, frag.sequ("right"), 1)

                                ss_left  = right_ss[0:idx_sta]
                                ss_right = right_ss[idx_sto:]
                                right_ss = ss_left + frag.sequ("right") + ss_right

            rnam_tmp = self.rnams[i]
            ss_rnam  = f'{rnam_tmp: <{rwid}}'
            sss = "%s\t%s\t%s" %(ss_rnam, left_ss, right_ss)
            f.write(f"{sss}\n")

        f.write("\n")
        f.close()

        return 0


if __name__ == "__main__":
    import __ParseArgs
    import SingleSamAlign
    import SingleAlign
    import PairSamAlign
    import FragmentAlign

    import pprint
    pp = pprint.PrettyPrinter(indent = 4)

    args = __ParseArgs.parseArgs()

    rec = {}
    reads = []

    with open("Test_Data_SingleAlign.sam", "r") as f:
        for line in f:
            tmp = line.split("\t")

            flag1 = tmp[0].endswith("C3") or tmp[0].endswith("C3#")
            flag2 = tmp[0].endswith("C6") or tmp[0].endswith("C6#")
            if not (flag1 or flag2):
                continue

            idx = tmp[0].find("_")
            nam = tmp[0][:idx]

            if not nam in rec:
                reads.append(nam)
                rec[ nam ] = []

            rec[ nam ].append(line[:-1])

    # print(rec)
    # pp.pprint(rec)

    # for read_tmp in reads:
    #     id = read_tmp
    #     print(">%s" %(id))

    #     pairAln = PairSamAlign.PairSamAlign(rec[read_tmp], args)
    #     frag_ref = pairAln.toFragment()
    #     pp.pprint(frag_ref)

    #     frag = FragmentAlign.FragmentAlign(frag_ref, args)

    #     fragCluster = FragmentCluster(frag)
    #     print(f"Read name: {fragCluster.rnams}")
    #     pp.pprint(fragCluster.frags[0].right)

    #     print()

    read_tmp = reads[0]
    pairAln = PairSamAlign.PairSamAlign(rec[read_tmp], args)
    frag_ref = pairAln.toFragment()
    frag1 = FragmentAlign.FragmentAlign(frag_ref, args)
    
    frag_clus = FragmentCluster()
    if frag_clus.isEmpty():
        print("The fragment cluster is empty.")

    frag_clus.addToCluster(frag1)
    # if not frag_clus.isEmpty():
    #     print("the fragment cluster is not empty.")

    read_tmp = reads[1]
    pairAln = PairSamAlign.PairSamAlign(rec[read_tmp], args)
    frag_ref = pairAln.toFragment()
    frag2 = FragmentAlign.FragmentAlign(frag_ref, args)

    if frag_clus.fragOverlap(frag2, "virus"):
        print("Overlap!")

    if frag_clus.mergeStatus(frag2):
        print("Can be merged.")

    if frag_clus.virus_rev() == 1:
        print("The virus part is reversed.")
    else:
        print("The virus part is not reversed.")

    if frag_clus.virus_circ() == 1:
        print("The virus part is circular.")
    else:
        print("The virus part is not circular.")

    # Human distance:
    dis = frag_clus.frag_distance(frag2)
    print(f"Human distance: {dis}")

    # Virus distance:
    dis = frag_clus.frag_distance(frag2, "virus")
    print(f"Virus distance: {dis}")

    frag_clus.addToCluster(frag2)
    print(frag_clus.rnams)

    pp.pprint(frag_clus.frags_dict)

    frag_sort = frag_clus.sortFrags()
    print(frag_sort.rnams)

    bpScore = frag_clus.breakpointScore()
    pp.pprint(bpScore)

    homolog = frag_clus.homologCall()
    pp.pprint(homolog)

    ref_merge = frag_clus.mergeClusterFrags()
    pp.pprint(ref_merge)

    frag_clus.printFragCluster("FragCluster_res.txt")

