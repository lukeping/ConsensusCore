
from pbcore.io import FastaReader, FastaWriter

import ConsensusCore as cc
import numpy as np
import numpy.random as nr

def simSmrtBellReads(ccRng, seqParams, insert, readLength):
    fwdTpl = insert
    revTpl = cc.ReverseComplement(insert)

    reads = []
    while sum([len(r) for r in reads]) < readLength:
        if len(reads) % 2 == 0:
            tpl = fwdTpl
        else:
            tpl = revTpl
        simRead = cc.SimulateRead(seqParams, tpl, ccRng)

        # acc = cc.Align(tpl, simRead).Accuracy()
        # print "Acc: ", acc

        reads.append(simRead)

    # Truncate last read
    requiredLastReadLength = readLength - sum([len(r) for r in reads[-1]])
    reads[-1] = reads[-1][:requiredLastReadLength]
    return reads

def sampleInsertFromChromosome(chromosome, insertSize, circular=False):
    #
    # Returns tStart, tEnd, strand, insert
    #
    chromosomeLen = len(chromosome)
    strand = nr.randint(0, 2)
    start  = nr.randint(0, chromosomeLen)
    tStart = start
    tEnd   = ((tStart + insertSize) % chromosomeLen) if circular \
             else min(chromosomeLen, tStart + insertSize)
    insert_ = chromosome[tStart:tEnd] if tStart <= tEnd \
              else genome[:tEnd] + genome[tStart:]
    insert = cc.ReverseComplement(insert_) if strand \
             else insert_

    return tStart, tEnd, strand, insert


def simExperiment(seed, seqParams,
                  meanReadLength,
                  meanInsertSize,
                  genome,
                  numSmrtBells,
                  circular=True):
    nr.seed(seed)
    ccRng = cc.RandomNumberGenerator(seed)

    for moleculeId in xrange(numSmrtBells):
        chromosomeId = nr.randint(0, len(genome))
        chromosome = genome[chromosomeId]

        # No model for readLengths, insert lengths as yet
        readLength = meanReadLength
        insertSize = meanInsertSize

        tStart, tEnd, strand, insert = sampleInsertFromChromosome(chromosome, insertSize, circular)
        if len(insert) < 500: continue
        smrtBellReads = simSmrtBellReads(ccRng, seqParams, insert, readLength)

        for subreadId, subread in enumerate(smrtBellReads):
            yield moleculeId, subreadId, (strand + subreadId) % 2, chromosomeId, tStart, tEnd, subread


if __name__ == '__main__':
    fw = FastaWriter("/tmp/out.fa")
    #genome = [r.sequence for r in FastaReader("~/Data/lambdaNEB.fa")]
    genome = [r.sequence for r in FastaReader("~/Data/Diploid/diploidLambda.fa")]
    seqParams = cc.SequencingParameters.C2()
    for item in simExperiment(42, seqParams, 4000, 1000, genome, 50000, False):
        fw.writeRecord("r@" + "_".join(map(str, item[:-1])),
                       item[-1])
