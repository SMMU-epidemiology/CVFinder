#!/usr/bin/env python

import PairSamAlign
import FragmentAlign
import FragmentCluster
import IdeogramAlign
import __Utils
import __ParseArgs


##### Auxilary functions:
#

readSam        = __Utils.readSam
getChromosomes = __Utils.getChromosomes

#
#####


if __name__ == "__main__":
    import pprint
    pp = pprint.PrettyPrinter(indent = 4)

    args = __ParseArgs.parseArgs()

    print("===== Collecting chimeric reads ... =====")
    
    rids, sam, virus_len = readSam(args.infile, args.virus)
    if virus_len != -1:
        args.vlen = virus_len
    
    print("=====            Done!              =====")
    
    if len(rids) == 0:
        print("No chimeric reads were detected!")
        exit(1)

    ideogram = IdeogramAlign.IdeogramAlign(args)

    for rid_tmp in rids:
        pair_align = PairSamAlign.PairSamAlign(
            sam[rid_tmp], 
            args
        )

        if pair_align.alignQC(args.mapq):
            continue

        frag_ref = pair_align.toFragment()
        frag = FragmentAlign.FragmentAlign(frag_ref, args)

        if frag.is_dummy():
            continue

        ideogram.addFragment(frag)

    chrs_ref    = getChromosomes(args)
    chrs        = chrs_ref["human"]
    orientation = ideogram.orientation

    chrs        = sorted(chrs)
    orientation = sorted(orientation)

    for orient in orientation:
        for chro in chrs:
            if ideogram.nfrags(orient, chro) < 1:
                continue

            clus_vec = ideogram.makeCluster(orient, chro)
            nclus    = len(clus_vec)

            if nclus == 0:
                continue

            clus_vector = []
            if args.refine == 1:
                for i in range(0, nclus):
                    clus_tmp = clus_vec[i]

                    vec  = clus_tmp.refineCluster()
                    nvec = len(vec)

                    for j in range(0, nvec):
                        vec[j].sortFrags()
                        clus_vector.append(vec[j])
            else:
                for i in range(0, nclus):
                    clus_vector.append(clus_vec[i])

            nclus = len(clus_vector)

            if nclus == 1:
                if clus_vector[0].nfrags < args.cluster_size:
                    continue

                clus_vector[0].printFragCluster(args.out)
                continue

            for i in range(0, nclus-1):
                clus_current = clus_vector[i]

                if clus_current.mergeStat() == 1:
                    continue

                for j in range(i+1, nclus-1):
                    clus_tmp = clus_vector[j]

                    if clus_tmp.mergeStat() == 1:
                        continue

                    if clus_current.mergeClusterStatus(clus_tmp) == 1:
                        clus_current = clus_current.mergeCluster(clus_tmp)
                        clus_tmp.mergeStat(1)

                        continue

                if clus_current.nfrags < args.cluster_size:
                    continue

                clus_current.sortFrags()
                clus_current.printFragCluster(args.out)

            if clus_vector[nclus-1].mergeStat() == 0:
                if clus_vector[nclus-1].nfrags < args.cluster_size:
                    continue

                clus_vector[nclus-1].printFragCluster(args.out)

