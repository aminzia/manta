#
# Manta
# Copyright (c) 2013-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

"""
Manta SV discovery workflow
"""


import os.path
import shutil
import sys

# add script path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(scriptDir))

# add pyflow path:
# TODO: get a more robust link to the pyflow dir at config time:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from pyflow import WorkflowRunner
from workflowUtil import checkFile, ensureDir, preJoin, which, \
                         getChromIntervals, getFastaChromOrderSize

from configureUtil import getIniSections,dumpIniSections



def getVersion() :
    return "@MANTA_VERSION@"


__version__ = getVersion()



def cleanId(input_id) :
    """
    filter id so that it's safe to use as a pyflow indentifier
    """
    import re
    return re.sub(r'([^a-zA-Z0-9_\-])', "_", input_id)



class GenomeSegment(object) :
    """
    organizes all variables which can change
    with each genomic segment.

    The genomic segment is defined by:

    1. chromosome
    2. begin position (1-indexed closed)
    3. end position (1-indexed closed)
    4. chromosome segment (ie. bin) number (0-indexed)
    """

    def __init__(self,chromIndex,chromLabel,beginPos,endPos,binId,genomeRegion) :
        """
        arguments are the 4 genomic interval descriptors detailed in class documentation
        """
        self.chromLabel = chromLabel
        self.bamRegion = chromLabel + ':' + str(beginPos) + '-' + str(endPos)
        self.binId = binId
        self.binStr = str(binId).zfill(4)
        self.id = chromLabel + "_" + self.binStr

        regionId=cleanId(chromLabel)
        if genomeRegion is not None :
            if genomeRegion['start'] is not None :
                regionId += "-"+str(genomeRegion['start'])
                if genomeRegion['end'] is not None :
                    regionId += "-"+str(genomeRegion['end'])
        self.pyflowId = "chromId_%s_%s_%s" % (str(chromIndex).zfill(3), regionId, self.binStr)



def getNextGenomeSegment(params) :
    """
    generator which iterates through all genomic segments and
    returns a segmentValues object for each one.
    """
    if params.genomeRegionList is None :
        for segval in getChromIntervals(params.chromOrder,params.chromSizes,params.scanSize) :
            yield GenomeSegment(*segval)
    else :
        for genomeRegion in params.genomeRegionList :
            for segval in getChromIntervals(params.chromOrder,params.chromSizes,params.scanSize, genomeRegion) :
                yield GenomeSegment(*segval)



def runStats(self,taskPrefix="",dependencies=None) :

    statsPath=self.paths.getStatsPath()

    cmd = [ self.params.mantaStatsBin ]
    cmd.extend(["--output-file",statsPath])
    for bamPath in self.params.normalBamList :
        cmd.extend(["--align-file",bamPath])
    for bamPath in self.params.tumorBamList :
        cmd.extend(["--tumor-align-file",bamPath])

    statsTask = self.addTask(preJoin(taskPrefix,"generateStats"),cmd,dependencies=dependencies)

    nextStepWait = set()
    nextStepWait.add(statsTask)

    # summarize stats for humans, no need for follow-up tasks to wait for this:
    cmd  = self.params.mantaStatsSummaryBin
    cmd += " --align-stats " + statsPath
    cmd += " > " + self.paths.getStatsSummaryPath()
    self.addTask(preJoin(taskPrefix,"summarizeStats"),cmd,dependencies=statsTask)

    return nextStepWait



def runDepth(self,taskPrefix="",dependencies=None) :
    """
    estimate chrom depth
    """

    bamFile=""
    if len(self.params.normalBamList) :
        bamFile = self.params.normalBamList[0]
    elif len(self.params.tumorBamList) :
        bamFile = self.params.tumorBamList[0]
    else :
        return set()


    cmd  = "%s -E %s" % (sys.executable, self.params.mantaChromDepth)
    cmd += " --bam '%s'" % (bamFile)
    cmd += " > %s" % (self.paths.getChromDepth())

    nextStepWait = set()
    nextStepWait.add(self.addTask(preJoin(taskPrefix,"estimateChromDepth"),cmd,dependencies=dependencies))

    return nextStepWait



def runLocusGraph(self,taskPrefix="",dependencies=None):
    """
    Create the full SV locus graph
    """

    statsPath=self.paths.getStatsPath()
    graphPath=self.paths.getGraphPath()
    graphStatsPath=self.paths.getGraphStatsPath()

    graphFilename=os.path.basename(graphPath)
    tmpGraphDir=os.path.join(self.params.workDir,graphFilename+".tmpdir")
    dirTask=self.addTask(preJoin(taskPrefix,"makeTmpDir"), "mkdir -p "+tmpGraphDir, dependencies=dependencies, isForceLocal=True)

    tmpGraphFiles = []
    graphTasks = set()

    for gseg in getNextGenomeSegment(self.params) :

        tmpGraphFiles.append(os.path.join(tmpGraphDir,graphFilename+"."+gseg.id+".bin"))
        graphCmd = [ self.params.mantaGraphBin ]
        graphCmd.extend(["--output-file", tmpGraphFiles[-1]])
        graphCmd.extend(["--align-stats",statsPath])
        graphCmd.extend(["--region",gseg.bamRegion])
        graphCmd.extend(["--min-candidate-sv-size", self.params.minCandidateVariantSize])
        graphCmd.extend(["--ref",self.params.referenceFasta])
        for bamPath in self.params.normalBamList :
            graphCmd.extend(["--align-file",bamPath])
        for bamPath in self.params.tumorBamList :
            graphCmd.extend(["--tumor-align-file",bamPath])

        if self.params.isHighDepthFilter :
            graphCmd.extend(["--chrom-depth", self.paths.getChromDepth()])

        if self.params.isIgnoreAnomProperPair :
            graphCmd.append("--ignore-anom-proper-pair")

        graphTaskLabel=preJoin(taskPrefix,"makeLocusGraph_"+gseg.pyflowId)
        graphTasks.add(self.addTask(graphTaskLabel,graphCmd,dependencies=dirTask,memMb=self.params.estimateMemMb))

    mergeCmd = [ self.params.mantaGraphMergeBin ]
    mergeCmd.extend(["--output-file", graphPath])
    for gfile in tmpGraphFiles :
        mergeCmd.extend(["--graph-file", gfile])

    mergeTask = self.addTask(preJoin(taskPrefix,"mergeLocusGraph"),mergeCmd,dependencies=graphTasks,memMb=self.params.mergeMemMb)

    # Run a separate process to rigorously check that the final graph is valid, the sv candidate generators will check as well, but
    # this makes the check much more clear:

    checkCmd = [ self.params.mantaGraphCheckBin ]
    checkCmd.extend(["--graph-file", graphPath])
    checkTask = self.addTask(preJoin(taskPrefix,"checkLocusGraph"),checkCmd,dependencies=mergeTask,memMb=self.params.mergeMemMb)

    rmGraphTmpCmd = "rm -rf " + tmpGraphDir
    rmTask=self.addTask(preJoin(taskPrefix,"rmGraphTmp"),rmGraphTmpCmd,dependencies=mergeTask)

    graphStatsCmd  = self.params.mantaGraphStatsBin
    graphStatsCmd += " --global"
    graphStatsCmd += " --graph-file " + graphPath
    graphStatsCmd += " >| " + graphStatsPath

    graphStatsTask = self.addTask(preJoin(taskPrefix,"locusGraphStats"),graphStatsCmd,dependencies=mergeTask,memMb=self.params.mergeMemMb)

    nextStepWait = set()
    nextStepWait.add(checkTask)
    return nextStepWait



def runHyGen(self, taskPrefix="", dependencies=None) :
    """
    Run hypothesis generation on each SV locus
    """

    import copy

    statsPath=self.paths.getStatsPath()
    graphPath=self.paths.getGraphPath()
    hygenDir=self.paths.getHyGenDir()

    dirTask=self.addTask(preJoin(taskPrefix,"makeHyGenDir"), "mkdir -p "+ hygenDir, dependencies=dependencies, isForceLocal=True)

    isSomatic = (len(self.params.normalBamList) and len(self.params.tumorBamList))

    hyGenMemMb = self.params.hyGenLocalMemMb
    if self.getRunMode() == "sge" :
        hyGenMemMb = self.params.hyGenSGEMemMb

    hygenTasks=set()
    candidateVcfPaths = []
    diploidVcfPaths = []
    somaticVcfPaths = []

    edgeLogPaths = []

    for binId in range(self.params.nonlocalWorkBins) :
        binStr = str(binId).zfill(4)
        candidateVcfPaths.append(self.paths.getHyGenCandidatePath(binStr))
        diploidVcfPaths.append(self.paths.getHyGenDiploidPath(binStr))
        if isSomatic :
            somaticVcfPaths.append(self.paths.getHyGenSomaticPath(binStr))

        edgeLogPaths.append(self.paths.getHyGenEdgeLogPath(binStr))

        hygenCmd = [ self.params.mantaHyGenBin ]
        hygenCmd.extend(["--align-stats",statsPath])
        hygenCmd.extend(["--graph-file",graphPath])
        hygenCmd.extend(["--bin-index", str(binId)])
        hygenCmd.extend(["--bin-count", str(self.params.nonlocalWorkBins)])
        hygenCmd.extend(["--min-candidate-sv-size", self.params.minCandidateVariantSize])
        hygenCmd.extend(["--min-scored-sv-size", self.params.minScoredVariantSize])
        hygenCmd.extend(["--ref",self.params.referenceFasta])
        hygenCmd.extend(["--candidate-output-file", candidateVcfPaths[-1]])
        hygenCmd.extend(["--diploid-output-file", diploidVcfPaths[-1]])
        hygenCmd.extend(["--min-qual-score", self.params.minDiploidVariantScore])
        hygenCmd.extend(["--min-pass-gt-score", self.params.minPassGTScore])
        if isSomatic :
            hygenCmd.extend(["--somatic-output-file", somaticVcfPaths[-1]])
            hygenCmd.extend(["--min-somatic-score", self.params.minSomaticScore])
            hygenCmd.extend(["--min-pass-somatic-score", self.params.minPassSomaticScore])

        if self.params.isHighDepthFilter :
            hygenCmd.extend(["--chrom-depth", self.paths.getChromDepth()])

        hygenCmd.extend(["--edge-runtime-log", edgeLogPaths[-1]])

        for bamPath in self.params.normalBamList :
            hygenCmd.extend(["--align-file", bamPath])
        for bamPath in self.params.tumorBamList :
            hygenCmd.extend(["--tumor-align-file", bamPath])
        
        if self.params.isIgnoreAnomProperPair :
            hygenCmd.append("--ignore-anom-proper-pair")
        if self.params.isRNA :
            hygenCmd.append("--rna")
            
        hygenTaskLabel=preJoin(taskPrefix,"generateCandidateSV_"+binStr)
        hygenTasks.add(self.addTask(hygenTaskLabel,hygenCmd,dependencies=dirTask, memMb=hyGenMemMb))

    nextStepWait = copy.deepcopy(hygenTasks)

    def getVcfSortCmd(vcfPaths, outPath) :
        cmd  = "%s -E %s -u " % (sys.executable,self.params.mantaSortVcf)
        cmd += " ".join(vcfPaths)
        cmd += " | %s -c > %s && %s -p vcf %s" % (self.params.bgzipBin, outPath, self.params.tabixBin, outPath)
        return cmd

    def sortVcfs(pathList, outPath, label) :
        if len(pathList) == 0 : return

        sortCmd = getVcfSortCmd(pathList,outPath)
        sortLabel=preJoin(taskPrefix,label)
        nextStepWait.add(self.addTask(sortLabel,sortCmd,dependencies=hygenTasks))

    sortVcfs(candidateVcfPaths, self.paths.getSortedCandidatePath(), "sortCandidateSV")
    sortVcfs(diploidVcfPaths, self.paths.getSortedDiploidPath(), "sortDiploidSV")
    sortVcfs(somaticVcfPaths, self.paths.getSortedSomaticPath(), "sortSomaticSV")

    # sort edge logs:
    edgeSortLabel=preJoin(taskPrefix,"sortEdgeLogs")
    edgeSortCmd="sort -rnk2 " + " ".join(edgeLogPaths) + " >| " + self.paths.getSortedEdgeLogPath()
    self.addTask(edgeSortLabel, edgeSortCmd, dependencies=hygenTasks, isForceLocal=True)

    return nextStepWait



class PathInfo:
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        self.params = params

    def getStatsPath(self) :
        return os.path.join(self.params.workDir,"alignmentStats.xml")

    def getStatsSummaryPath(self) :
        return os.path.join(self.params.statsDir,"alignmentStatsSummary.txt")

    def getChromDepth(self) :
        return os.path.join(self.params.workDir,"chromDepth.txt")

    def getGraphPath(self) :
        return os.path.join(self.params.workDir,"svLocusGraph.bin")

    def getHyGenDir(self) :
        return os.path.join(self.params.workDir,"svHyGen")

    def getHyGenCandidatePath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"candidateSV.%s.vcf" % (binStr))

    def getSortedCandidatePath(self) :
        return os.path.join(self.params.variantsDir,"candidateSV.vcf.gz")

    def getHyGenDiploidPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"diploidSV.%s.vcf" % (binStr))

    def getSortedDiploidPath(self) :
        return os.path.join(self.params.variantsDir,"diploidSV.vcf.gz")

    def getHyGenSomaticPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"somaticSV.%s.vcf" % (binStr))

    def getSortedSomaticPath(self) :
        return os.path.join(self.params.variantsDir,"somaticSV.vcf.gz")

    def getHyGenEdgeLogPath(self, binStr) :
        return os.path.join(self.getHyGenDir(),"edgeRuntimeLog.%s.txt" % (binStr))

    def getSortedEdgeLogPath(self) :
        return os.path.join(self.params.workDir,"edgeRuntimeLog.txt")

    def getGraphStatsPath(self) :
        return os.path.join(self.params.statsDir,"svLocusGraphStats.tsv")



class MantaWorkflow(WorkflowRunner) :
    """
    Manta SV discovery workflow
    """

    def __init__(self,params,iniSections) :

        # clear out some potentially destabilizing env variables:
        clearList = [ "PYTHONPATH", "PYTHONHOME"]
        for key in clearList :
            if key in os.environ :
                del os.environ[key]

        self.params=params
        self.iniSections=iniSections

        # format bam lists:
        if self.params.normalBamList is None : self.params.normalBamList = []
        if self.params.tumorBamList is None : self.params.tumorBamList = []

        # make sure run directory is setup:
        self.params.runDir=os.path.abspath(self.params.runDir)
        ensureDir(self.params.runDir)

        # everything that's not intended to be a final result should dump directories/files in workDir
        self.params.workDir=os.path.join(self.params.runDir,"workspace")
        ensureDir(self.params.workDir)

        # all finalized pretty results get transfered to resultsDir
        self.params.resultsDir=os.path.join(self.params.runDir,"results")
        ensureDir(self.params.resultsDir)
        self.params.statsDir=os.path.join(self.params.resultsDir,"stats")
        ensureDir(self.params.statsDir)
        self.params.variantsDir=os.path.join(self.params.resultsDir,"variants")
        ensureDir(self.params.variantsDir)
#         self.params.reportsDir=os.path.join(self.params.resultsDir,"reports")
#         ensureDir(self.params.reportsDir)

        indexRefFasta=self.params.referenceFasta+".fai"

        if self.params.referenceFasta is None:
            raise Exception("No reference fasta defined.")
        else:
            checkFile(self.params.referenceFasta,"reference fasta")
            checkFile(indexRefFasta,"reference fasta index")

        # read fasta index
        (self.params.chromOrder,self.params.chromSizes) = getFastaChromOrderSize(indexRefFasta)

        # sanity check some parameter typing:
        MEGABASE = 1000000
        self.params.scanSize = int(self.params.scanSizeMb) * MEGABASE
        self.params.nonlocalWorkBins = int(self.params.nonlocalWorkBins)

        self.paths = PathInfo(self.params)

        self.params.isHighDepthFilter = (not (self.params.isExome or self.params.isRNA))
        self.params.isIgnoreAnomProperPair = (self.params.isRNA)



    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Manta workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Manta workflow version: %s" % (__version__))

        graphTaskDependencies = set()

        if not self.params.useExistingAlignStats :
            statsTasks = runStats(self)
            graphTaskDependencies |= statsTasks

        if not ((not self.params.isHighDepthFilter) or self.params.useExistingChromDepths) :
            depthTasks = runDepth(self)
            graphTaskDependencies |= depthTasks

        graphTasks = runLocusGraph(self,dependencies=graphTaskDependencies)

        hygenTasks = runHyGen(self,dependencies=graphTasks)
