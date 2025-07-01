#!/usr/bin/env python

import cigar
import re
import math
import subprocess


def intersect(lis1, lis2):
    res = [v for v in lis1 if v in lis2]
    return res


def distance(args, spec = None):
    if spec == None:
        spec = "human"
    
    if spec == "human":
        return args.distance
    else:
        return args.virus_thre


# def compSeq(sequence):
#     comp_dict = {
#         "A":"T",
#         "T":"A",
#         "G":"C",
#         "C":"G",
#         "a":"t",
#         "t":"a",
#         "g":"c",
#         "c":"g",
#         "-":"-",
#         "N":"N",
#     }

#     sequence_list = list(sequence)
#     sequence_list = [comp_dict[base] for base in sequence_list]
#     string = ''.join(sequence_list)

#     return string


def compSeq(sequ):
    sequ = sequ.upper()
    
    sequ = sequ.replace("A", "1")
    sequ = sequ.replace("T", "A")
    sequ = sequ.replace("1", "T")

    sequ = sequ.replace("C", "1")
    sequ = sequ.replace("G", "C")
    sequ = sequ.replace("1", "G")

    return sequ


def revCompSequ(sequ):
    sequ = sequ[::-1]
    sequ = compSeq(sequ)
    return sequ


def bioExplainCigar(cigar_str):
    ref = cigar.Cigar(cigar_str)
    res = {}

    for i, ref_tmp in enumerate(ref.items()):
        res[i] = {
            "nbase" : ref_tmp[0],
            "symbol": ref_tmp[1],
        }

    return res


def reverseCigar(cigar):
    cigar_ref = bioExplainCigar(cigar)

    keys = list(cigar_ref.keys())
    keys_rev = keys[::-1]

    res = {}
    for key_tmp, rkey_tmp in zip(keys, keys_rev):
        res[key_tmp] = cigar_ref[rkey_tmp]

    return res


def set_reverse_cigar(cigar):
    cigar_rev = reverseCigar(cigar)
    keys = cigar_rev.keys()

    res = ""
    for key_tmp in keys:
        nbase  = cigar_rev[key_tmp]["nbase"]
        symbol = cigar_rev[key_tmp]["symbol"]
        res += str(nbase)
        res += symbol

    return res


def pruneSequ(sequ, qual, cigar):
    cigar_ref = bioExplainCigar(cigar)
    keys = cigar_ref.keys()

    idx = 0
    res = {}
    for key_tmp in keys:
        pat = re.compile("S")

        if re.search(pat, cigar_ref[key_tmp]["symbol"]) != None:
            idx += cigar_ref[key_tmp]["nbase"]

        if cigar_ref[key_tmp]["symbol"] == "M":
            nbase = cigar_ref[key_tmp]["nbase"]

            res["sequ"] = sequ[idx:idx+nbase]
            res["qual"] = qual[idx:idx+nbase]
            break

    return res


"""
ref = dict(
    start = 0,
    stop  = 0,
    sequ  = "NA",
    qual  = "NA",
)
"""

def is_dummy(ref):
    flag = 0

    if ref["sequ"] == "NA" and ref["qual"] == "NA":
        flag = 1

    return flag


"""
seg = dict(
    chro   = "NA",
    strand = "NA",
    start  = 0,
    stop   = 0,
    sequ   = "NA",
    qual   = "NA",
    rev    = 0,
    circ   = 0,
)
"""

def segMap(seg, vlen = 3215):
    res = {}

    sequ = seg["sequ"]
    qual = seg["qual"]

    sta  = seg["start"]
    sto  = seg["stop"]
    len  = None

    if seg["circ"] == 0:
        len = abs(sto - sta) + 1
    else:
        if seg["rev"] == 0:
            len = (vlen - sta + 1) + sto
        else:
            len = (vlen - sto + 1) + sta

    #   my $rev_stat = $seg->{rev};
    #   $rev_stat = 1 if $sta > $sto;

    cods = []
    if seg["circ"] == 0:
        if seg["rev"] == 0:
            cods = list(range(sta, sto+1))
        else:
            cods = list(range(sta, sto-1, -1))
    else:
        if seg["rev"] == 0:
            cods = list(range(sta, vlen+1))
            cods.extend(range(1, sto+1))
        else:
            cods = list(range(sto, vlen+1))
            cods.extend(range(1, sta+1))
            cods = cods[::-1]

    for i in range(0, len):
        res[ cods[i] ] = {}
        res[ cods[i] ]["let"] = sequ[i:i+1]
        res[ cods[i] ]["qua"] = qual[i:i+1]

    return res


def mergeSegs(seg1, seg2, vlen = 3215, circular_stat = 0):
    res = {
        "chro"   : seg1["chro"],
        "strand" : seg1["strand"],
        "start"  : 0,
        "stop"   : 0,
        "sequ"   : "NA",
        "qual"   : "NA",
        "rev"    : seg1["rev"],
        "circ"   : 0,
    }

    gap = {
       # 0: does not have gap; 1: have gap; 
       "gap_flag"   : 0,
       "gap_chro"   : "NA",
       "gap_strand" : ".",
       "gap_start"  : 0,
       "gap_stop"   : 0,
       "gap_rev"    : 0,
    }

    if is_dummy(seg2) == 1:
        return seg1

    if is_dummy(seg1) == 1:
        return seg2

    if seg1["circ"] == 1 or seg2["circ"] == 1:
        res["circ"] = 1

    sta1 = seg1["start"]
    sto1 = seg1["stop"]
    sta2 = seg2["start"]
    sto2 = seg2["stop"]

    cods = []
    sta = sta1
    sto = sto1
    rev_stat = seg1["rev"]
    if rev_stat == 0:
        if circular_stat == 0:
            if sta2 < sta:
                sta = sta2

            if sto2 > sto:
                sto = sto2

            cods = list(range(sta, sto+1))
        else:
            if seg1["circ"] == 0 and seg2["circ"] == 0:
                cod1 = (sta1 + sto1)/2
                cod2 = (sta2 + sto2)/2

                flag1 = 0
                flag2 = 0

                if vlen - cod1 < cod1 - 1:
                    flag1 = 1

                if vlen - cod2 < cod2 - 1:
                    flag2 = 1

                if flag1 == 0 and flag2 == 1:
                    if abs(cod1 - cod2) < (vlen - cod2 + 1) + cod1:
                        if sta2 < sta:
                            sta = sta2

                        if sto2 > sto:
                            sto = sto2

                        cods = list(range(sta, sto+1))
                    else:
                        res["circ"] = 1
                        sta = sta2
                        sto = sto1

                        cods = list(range(sta, vlen+1))
                        cods.extend(range(1, sto+1))
                elif flag1 == 1 and flag2 == 0:
                    if abs(cod1 - cod2) < (vlen - cod1 + 1) + cod2:
                        if sta2 < sta:
                            sta = sta2

                        if sto2 > sto:
                            sto = sto2

                        cods = list(range(sta, sto+1))
                    else:
                        res["circ"] = 1
                        sta = sta1
                        sto = sto2

                        cods = list(range(sta, vlen+1))
                        cods.extend(range(1, sto+1))
                else:
                    if sta2 < sta:
                        sta = sta2

                    if sto2 > sto:
                        sto = sto2

                    cods = list(range(sta, sto+1))
            elif seg1["circ"] == 1 and seg2["circ"] == 0:
                if vlen - sta2 < sta2 - 1:
                    if sta2 < sta:
                        sta = sta2

                    sto = sto1
                else:
                    sta = sta1

                    if sto2 > sto:
                        sto = sto2

                cods = list(range(sta, vlen+1))
                cods.extend(range(1, sto+1))
            elif seg1["circ"] == 0 and seg2["circ"] == 1:
                # Unnecessary
                res["circ"] = 1

                if vlen - sta < sta - 1:
                    if sta2 < sta:
                        sta = sta2

                    sto = sto2
                else:
                    sta = sta2

                    if sto2 > sto1:
                        sto = sto2

                cods = list(range(sta, vlen+1))
                cods.extend(range(1, sto+1))
            else:
                if sta2 < sta:
                    sta = sta2

                if sto2 > sto:
                    sto = sto2

                cods = list(range(sta, vlen+1))
                cods.extend(range(1, sto+1))
    else:
        if circular_stat == 0:
            if sta2 > sta:
                sta = sta2

            if sto2 < sto:
                sto = sto2

            cods = list(range(sto, sta+1))[::-1]
        else:
            if seg1["circ"] == 0 and seg2["circ"] == 0:
                cod1 = (sta1 + sto1)/2
                cod2 = (sta2 + sto2)/2

                flag1 = 0
                flag2 = 0

                if vlen - cod1 < cod1 - 1:
                    flag1 = 1

                if vlen - cod2 < cod2 - 1:
                    flag2 = 1

                if flag1 == 0 and flag2 == 1:
                    if abs(cod1 - cod2) < (vlen - cod2 + 1) + cod1:
                        if sta2 > sta:
                            sta = sta2

                        if sto2 < sto:
                            sto = sto2

                        cods = list(range(sto, sta+1))[::-1]
                    else:
                        res["circ"] = 1
                        sta = sta1
                        sto = sto2

                        cods = list(range(sto, vlen+1))
                        cods.extend(range(1, sta+1))
                        cods = cods[::-1]
                elif flag1 == 1 and flag2 == 0:
                    if abs(cod1 - cod2) < (vlen - cod1 + 1) + cod2:
                        if sta2 > sta:
                            sta = sta2

                        if sto2 < sto:
                            sto = sto2

                        cods = list(range(sto, sta+1))[::-1]
                    else:
                        res["circ"] = 1
                        sta = sta2
                        sto = sto1

                        cods = list(range(sto, vlen+1))
                        cods.extend(range(1, sta+1))
                        cods = cods[::-1]
                else:
                    if sta2 > sta:
                        sta = sta2

                    if sto2 < sto:
                        sto = sto2

                    cods = list(range(sto, sta+1))[::-1]
            elif seg1["circ"] == 1 and seg2["circ"] == 0:
                if vlen - sta2 < sta2 - 1:
                    sta = sta1

                    if sto2 < sto:
                        sto = sto2
                else:
                    if sta2 > sta:
                        sta = sta2

                    sto = sto1

                cods = list(range(sto, vlen+1))
                cods.extend(range(1, sta+1))
                cods = cods[::-1]
            elif seg1["circ"] == 0 and seg2["circ"] == 1:
                # Unnecessary
                res["circ"] = 1

                if vlen - sta < sta - 1:
                    sta = sta2

                    if sto2 < sto:
                        sto = sto2
                else:
                    if sta2 > sta:
                        sta = sta2

                    sto = sto2

                cods = list(range(sto, vlen+1))
                cods.extend(range(1, sta+1))
                cods = cods[::-1]
            else:
                if sta2 > sta:
                    sta = sta2

                if sto2 < sto:
                    sto = sto2

                cods = list(range(sto, vlen+1))
                cods.extend(range(1, sta+1))
                cods = cods[::-1]

    res["start"] = cods[0]
    res["stop"]  = cods[-1]

    seg_ref1 = segMap(seg1, vlen)
    seg_ref2 = segMap(seg2, vlen)

    sequ = ""
    qual = ""
    cod_buf = []

    for cod_tmp in cods:
        let = None
        qua = None

        if (not cod_tmp in seg_ref1) and (not cod_tmp in seg_ref2):
            let = "-"
            qua = "-"

            gap["gap_flag"]   = 1
            gap["gap_chro"]   = seg1["chro"]
            gap["gap_strand"] = seg1["strand"]
            gap["gap_rev"]    = seg1["rev"]

            cod_buf.append(cod_tmp)
        elif (cod_tmp in seg_ref1) and (not cod_tmp in seg_ref2): 
            if seg_ref1[cod_tmp] == "-":
                let = "-"
                qua = "-"
            else:
                let = seg_ref1[cod_tmp]["let"]
                qua = seg_ref1[cod_tmp]["qua"]
        elif (not cod_tmp in seg_ref1) and (cod_tmp in seg_ref2):
            if seg_ref2[cod_tmp] == "-":
                let = "-"
                qua = "-"
            else:
                let = seg_ref2[cod_tmp]["let"]
                qua = seg_ref2[cod_tmp]["qua"]
        else:
            let1 = seg_ref1[cod_tmp]["let"]
            qua1 = seg_ref1[cod_tmp]["qua"]
            let2 = seg_ref2[cod_tmp]["let"]
            qua2 = seg_ref2[cod_tmp]["qua"]

            if let1 == "-" and let2 == "-":
                let = "-"
                qua = "-"
            elif let1 == "-" and let2 != "-":
                let = let2
                qua = qua2
            elif let1 != "-" and let2 == "-":
                let = let1
                qua = qua1
            else:
                if qua1 >= qua2:
                    let = let1
                    qua = qua1
                else:
                    let = let2
                    qua = qua2

        sequ += let
        qual += qua

    res["sequ"] = sequ
    res["qual"] = qual

    if gap["gap_flag"] == 1:
        gap["gap_start"] = cod_buf[0]
        gap["gap_stop"]  = cod_buf[-1]

    return res, gap


def distanceVirus(cod1, cod2, vlen, circ_stat):
    if cod1 > cod2:
        cod1, cod2 = cod2, cod1

    if circ_stat == 0:
        return cod2 - cod1
    else:
        if cod2 - cod1 < vlen - cod2 + cod1:
            return cod2 - cod1
        else:
            return vlen - cod2 + cod1


# linear or circular
def distanceVirusPattern(cod1, cod2, vlen):
    if cod1 > cod2:
        cod1, cod2 = cod2, cod1

    if cod2 - cod1 < vlen - cod2 + cod1:
        return "linear"
    else:
        return "circular"


#############  Sequence Complexity  #############

WINDOWSIZE = 64
WINDOWSTEP = 32
WORDSIZE = 3
WINDOWSIZEARRAY = range(0, 62)
LOG62 = math.log(62)
ONEOVERLOG62 = 1/math.log(62)
POINTFIVE = 1/2

#################################################

def getArrayMean(lis = []):
    if len(lis) == 0:
        return 0
    else:
        return sum(lis)/len(lis)


# Dust Score:
def dustScore(seqn):
    res = None

    rest  = None
    steps = None
    vals  = []
    str   = None
    num   = None
    bynum = None

    length = len(seqn)
    if length <= WINDOWSIZE:
        rest  = length
        steps = 0
    else:
        steps = int((length - WINDOWSIZE) / WINDOWSTEP) + 1
        rest  = length - steps * WINDOWSTEP

        if not rest > WINDOWSTEP:
            rest  += WINDOWSTEP
            steps -= 1
    
    num = WINDOWSIZE - 2
    bynum = 1/num
    num -= 1
    mean = 0

    dustscore = None
    for i in range(0, steps):
        idx_sta = i*WINDOWSTEP
        idx_sto = idx_sta + WINDOWSIZE
        str = seqn[idx_sta:idx_sto]
        counts = dict()

        for i in WINDOWSIZEARRAY:
            ckey = str[i:i+3]

            if not ckey in counts:
                counts[ckey] = 0
            
            counts[ckey] += 1
        
        dustscore = 0
        for val in counts.values():
            dustscore += val * (val - 1) * POINTFIVE
        
        vals.append(dustscore * bynum)
    
    # last step
    if rest > 5:
        idx_sta = steps * WINDOWSTEP
        idx_sto = idx_sta + rest
        str = seqn[idx_sta:idx_sto]

        counts = {}
        num = rest - 2

        for i in range(0,  num):
            ckey = str[i:i+3]

            if not ckey in counts:
                counts[ckey] = 0
            
            counts[ckey] += 1
        
        dustscore = 0
        for val in counts.values():
            dustscore += (val * (val - 1) * POINTFIVE)
        
        vals.append((dustscore/(num-1)) * ((WINDOWSIZE-2)/num))
    else:
        #to assign a maximum score based on the scaling factor 100/31
        vals.append(31)
    
    res = getArrayMean(vals)
    res = int(res*100/31)

    return res


def entropyScore(seqn):
    res = None

    rest  = None
    steps = None
    vals  = []
    str   = None
    num   = None
    bynum = None

    length = len(seqn)
    if length <= WINDOWSIZE:
        rest  = length
        steps = 0
    else:
        steps = int((length - WINDOWSIZE)/WINDOWSTEP) + 1
        rest  = length - steps*WINDOWSTEP

        if not rest > WINDOWSTEP:
            rest += WINDOWSTEP
            steps -= 1
    
    num   = WINDOWSIZE-2
    bynum = 1/num
    num  -= 1
    mean  = 0

    entropyval = None
    for i in range(0, steps):
        idx_sta = i*WINDOWSTEP
        idx_sto = idx_sta + WINDOWSIZE
        str = seqn[idx_sta:idx_sto]

        counts = {}
        for i in WINDOWSIZEARRAY:
            ckey = str[i:i+3]

            if not ckey in counts:
                counts[ckey] = 0
            
            counts[ckey] += 1
        
        entropyval = 0
        for val in counts.values():
            entropyval -= (val*bynum)*math.log(val*bynum)
        
        vals.append(entropyval * ONEOVERLOG62)
    
    # last step
    if rest > 5:
        idx_sta = steps*WINDOWSTEP
        idx_sto = idx_sta + rest
        str = seqn[idx_sta:idx_sto]

        counts = {}
        num    = rest - 2

        for i in range(0, num):
            ckey = str[i:i+3]

            if not ckey in counts:
                counts[ckey] = 0
            
            counts[ckey] += 1
        
        entropyval = 0
        bynum = 1/num
        for val in counts.values():
            entropyval -= (val*bynum)*math.log(val*bynum)
        
        vals.append(entropyval/math.log(num))
    else:
        vals.append(0)
    
    res = getArrayMean(vals)
    res = round(res, 2)

    return res


# Retrieving all chromosome names from .sam file
# and taking last one as virus by default.
# Take name of .sam file as input
def getChromosomes(para):
    res = dict(
        human = [],
        virus = {},
    )

    filename = para.infile
    vnam     = para.virus

    if filename.endswith("sam"):
        p = open(filename, "r")

    if filename.endswith("bam"):
        p = subprocess.Popen(
            f"samtools view -h {filename}",
            shell  = True,
            stdout = subprocess.PIPE
        )

    vec = []
    while(True):
        if filename.endswith("sam"):
            line = p.readline()

        if filename.endswith("bam"):
            line = p.stdout.readline()
            line = line.decode()

        if line == "":
            break

        if not line.startswith("@"):
            break

        if not line.startswith("@SQ"):
            continue

        vec.append(line[:-1])

    if filename.endswith("sam"):
        p.close()

    if filename.endswith("bam"):
        p.terminate()


    for ss in vec:
        tmp  = ss.split("\t")
        chro = tmp[1][3:]
        clen = tmp[2][3:]

        if not chro.startswith("chr"):
            chro = "chr" + chro

        if chro == vnam:
            res["virus"] = dict(
                chro = chro,
                clen = clen,
            )

            continue

        res["human"].append(chro)

    return res


def isChimeric(vec, virus):
    flag_v = 0
    flag_h = 0

    for ss in vec:
        tmp = ss.split("\t")

        if tmp[2] == virus:
            flag_v = 1
        else:
            flag_h = 1

    if flag_v == 1 and flag_h == 1:
        return 1
    else:
        return 0


def readSam(filename, vnam):
    res  = {}
    rids = []

    rec = []
    rid_current = ""
    rid = ""

    pat = re.compile(vnam)
    virus_len = -1

    if filename.endswith("sam"):
        p = open(filename, "r")

    if filename.endswith("bam"):
        p = subprocess.Popen(
                f"samtools view -h {filename}",
                shell  = True,
                stdout = subprocess.PIPE
            )

    while(True):
        if filename.endswith("sam"):
            line = p.readline()

        if filename.endswith("bam"):
            line = p.stdout.readline()
            line = line.decode()

        if line == "":
            break

        if line.startswith("#"):
            continue

        if line.startswith("__END__"):
            break

        if line.startswith("@"):
            if re.search(pat, line) == None:
                continue
            else:
                tmp = line.split("\t")
                virus_len = int(tmp[2][3:-1])

                continue

        idx = line.find("\t")
        rid = line[:idx]

        if rid_current == "":
            rid_current = rid
            rec.append(line[:-1])

            continue

        if rid != rid_current:
            # Check if chimeric read or not
            flag = isChimeric(rec, vnam)

            if flag != 0:
                res[ rid_current ] = rec
                rids.append(rid_current)

            rid_current = rid
            rec = []
            rec.append(line[:-1])
            continue
        else:
            rec.append(line[:-1])
            continue

    if filename.endswith("sam"):
        p.close()

    if filename.endswith("bam"):
        p.terminate()

    flag = isChimeric(rec, vnam)

    if flag != 0:
        res[ rid ] = rec
        rids.append(rid)

    return rids, res, virus_len


def readInfo(filename):
    rec = []

    with open(filename, "r") as f:
        for line in f:
            rec.append(line[:-1])

    return rec


def readCluster(filename):
    res  = []
    clus = []

    with open(filename, "r") as f:
        line = f.readline()
        clus.append(line[:-1])

        for line in f:
            if line[:-1] == "":
                continue

            if line.startswith("Left"):
                res.append(clus)
                clus.append(line[:-1])

                continue

            clus.append(line[:-1])

    return res

