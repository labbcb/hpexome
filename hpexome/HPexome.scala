/*
* HPexome
* Copyright (C) 2017-2019 Lucas Lopes Cendes, Welliton de Souza
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._

class HPexome extends QScript {

  @Input(doc = "Reference Genome File", shortName = "R")
  var referenceFile: File = _

  @Input(doc = "BAM Files", shortName = "I")
  var bamFiles: Seq[File] = Nil

  @Input(doc = "Known Sites Files", shortName = "knownSites", required = true)
  var knownSitesFiles: List[File] = Nil

  @Input(doc = "Known Indels Files", shortName = "known", required = false)
  var knownIndelsFiles: List[File] = Nil

  @Input(doc = "dbSNP File", shortName = "dbsnp", required = true)
  var dbSnpFile: File = _

  @Input(doc = "Genomic intervals over which to operate", shortName = "L", required = false)
  var genomicIntervalsFile: Seq[File] = Nil

  @Argument(doc = "Number of data threads", shortName = "nt", required = false)
  var ntp: Int = _

  @Argument(doc = "Number of cpu threads", shortName = "nct", required = false)
  var nctp: Int = _

  @Argument(doc = "Scatter Count", shortName = "scattercount", required = false)
  var scatterCount: Int = _

  @Argument(doc = "name of the output folder", shortName = "outdir", required = false)
  var outputDir: String = _

  @Argument(doc = "A single output for several inputs in the haplotype caller", shortName = "unifiedVCF", required = false)
  var unifiedVcf: Boolean = false

  @Argument(doc = "Standard min confidence threshold for calling", shortName = "stand_call_conf", required = false)
  var standCallConf: Int = _

  @Argument(doc = "Minimum support to not prune paths in graph", shortName = "minPruning", required = false)
  var minPruning: Int = _

  @Argument(doc = "Output file name for unified VCF (--unified_vcf)", shortName = "filename", required = false)
  var outputFileName: String = "unified.vcf"

  @Argument(doc = "Maximum Java memory in gigabytes", shortName = "Xmx", required = false)
  var javaMemoryLimit: Double = 2

  @Output
  var unifiedVcfFile: File = outputFileName

  if (standCallConf == 0) {
    standCallConf = 30
  }
  if (minPruning == 0) {
    minPruning = 2
  }

  trait GenomeReference extends CommandLineGATK {
    this.reference_sequence = referenceFile
  }

  trait ParallelNCT extends CommandLineGATK {
    if (nctp > 1 && scatterCount <= 1) {
      this.nct = nctp
    }
  }

  trait GenomicIntervals extends CommandLineGATK {
    if (genomicIntervalsFile != null) {
      this.intervals = genomicIntervalsFile
    }
  }

  trait HCConfig extends HaplotypeCaller {
    this.stand_call_conf = standCallConf
    this.minPruning = minPruning
    this.dbsnp = dbSnpFile
  }

  def script() {
    if (outputDir == null) {
      outputDir = "."
    }

    var printReadsOutputFiles: Seq[File] = Nil

    for (bamFile <- bamFiles) {
      val realignerTargetCreator = new RealignerTargetCreator with GenomeReference
      realignerTargetCreator.input_file :+= bamFile
      if (knownIndelsFiles != null) {
        realignerTargetCreator.known = knownIndelsFiles
      }
      realignerTargetCreator.out = swapExt(outputDir, bamFile, "bam", "intervals")
      if (ntp > 1 && scatterCount <= 1) {
        realignerTargetCreator.nt = ntp
      }
      if (scatterCount > 1) {
        realignerTargetCreator.scatterCount = scatterCount
      }
      add(realignerTargetCreator)

      val indelRealigner = new IndelRealigner with GenomeReference with GenomicIntervals
      indelRealigner.input_file :+= bamFile
      if (knownIndelsFiles != null) {
        indelRealigner.known = knownIndelsFiles
      }
      indelRealigner.targetIntervals = realignerTargetCreator.out
      indelRealigner.out = swapExt(outputDir, bamFile, "bam", "realn.bam")
      if (scatterCount > 1) {
        indelRealigner.scatterCount = scatterCount
      }
      add(indelRealigner)

      val baseRecalibrator = new BaseRecalibrator with GenomeReference with ParallelNCT
      baseRecalibrator.input_file :+= indelRealigner.out
      if (knownSitesFiles != null) {
        baseRecalibrator.knownSites = knownSitesFiles
      }
      if (scatterCount > 1) {
        baseRecalibrator.scatterCount = scatterCount
      }
      baseRecalibrator.out = swapExt(outputDir, bamFile, "bam", "recal.cvs")
      add(baseRecalibrator)

      val printReads = new PrintReads with GenomeReference with ParallelNCT
      printReads.input_file :+= indelRealigner.out
      printReads.BQSR = baseRecalibrator.out
      printReads.out = swapExt(outputDir, bamFile, "bam", "recal.bam")
      if (scatterCount > 1) {
        printReads.scatterCount = scatterCount
      }
      add(printReads)

      if (!unifiedVcf) {
        val haplotypeCaller = new HaplotypeCaller with GenomeReference with HCConfig
        haplotypeCaller.input_file :+= printReads.out
        haplotypeCaller.out = swapExt(outputDir, bamFile, "bam", "HC.raw.vcf")
        if (scatterCount > 1) {
          haplotypeCaller.scatterCount = scatterCount
        }
        add(haplotypeCaller)
      }
      printReadsOutputFiles +:= printReads.out
    }

    if (unifiedVcf) {
      val haplotypeCaller = new HaplotypeCaller with GenomeReference with HCConfig with GenomicIntervals
      haplotypeCaller.input_file = printReadsOutputFiles
      haplotypeCaller.out = outputDir + "/" + unifiedVcfFile
      if (scatterCount > 1) {
        haplotypeCaller.scatterCount = scatterCount
      }
      haplotypeCaller.javaMemoryLimit = javaMemoryLimit
      add(haplotypeCaller)
    }
  }
}
