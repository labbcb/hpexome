/*
* hpexome
* Copyright (C) 2017-2019 Lucas Lopes Cendes
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
package script

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.ListWriterFunction
import java.io.File
import org.broadinstitute.gatk.queue.extensions.picard.PicardBamFunction

class Hpexome extends QScript {

  @Input(doc = "Reference file", shortName = "R")
  var reference: File = _

  @Input(doc = "One or more bamfiles", shortName = "I")
  var bam: Seq[File] = Nil

  @Input(doc = "Known Sites", shortName = "knownSites", required = true)
  var sites: List[File] = Nil

  @Input(doc = "Known indels", shortName = "known", required = false)
  var know: List[File] = Nil

  @Input(doc = "dbSNP file", shortName = "dbsnp", required = true)
  var dbsnip: File = _

  @Input(doc = "One or more genomic intervals over which to operate", shortName = "L", required = false)
  var intervals: File = _

  @Argument(doc = "Number of data threads", shortName = "nt", required = false)
  var ntp: Int = _

  @Argument(doc = "Number of cpu threads", shortName = "nct", required = false)
  var nctp: Int = _

  @Argument(doc = "Scatter Count", shortName = "scattercount", required = false)
  var scatter: Int = _

  @Argument(doc = "name of the output folder", shortName = "outdir", required = false)
  var outdir: String = _

  @Argument(doc = "A single output for several inputs in the haplotype caller", shortName = "unifiedVCF", required = false)
  var unifiedvcf: Boolean = false

  @Argument(doc = "Standard min confidence threshold for calling", shortName = "stand_call_conf", required = false)
  var standcallconf: Int = _

  @Argument(doc = "Minimum support to not prune paths in graph", shortName = "minPruning", required = false)
  var minpruning: Int = _

  @Argument(doc = "Output file name for unified VCF (--unified_vcf)", shortName = "o", required = false)
  var outputFileName: String = "unified.vcf"

  @Output
  var hapcall: File = outputFileName

  if (standcallconf == 0) {
    standcallconf = 30
  }
  if (minpruning == 0) {
    minpruning = 2
  }

  trait referenceseq extends CommandLineGATK {
    this.reference_sequence = reference
  }

  trait nctparallel extends CommandLineGATK {
    if (nctp > 1 && scatter <= 1) {
      this.nct = nctp
    }
  }

  trait intervals extends CommandLineGATK {
    if (intervals != null) {
      this.intervals = intervals
    }
  }

  trait hapcall extends HaplotypeCaller {
    this.stand_call_conf = standcallconf
    this.minPruning = minpruning
    this.dbsnp = dbsnip
  }

  def script() {
    if (outdir == null) {
      outdir = "."
    } else {
      Runtime.getRuntime().exec("mkdir " + outdir)
    }

    var inputhap: Seq[File] = Nil

    for (a <- bam) {
      val trealigner = new RealignerTargetCreator with referenceseq
      val indelrealigner = new IndelRealigner with referenceseq with intervals
      val brecalibrator = new BaseRecalibrator with referenceseq with nctparallel
      val printreads = new PrintReads with referenceseq with nctparallel

      trealigner.input_file :+= a
      if (know != null) {
        trealigner.known = know
      }
      trealigner.out = swapExt(outdir, a, "bam", "intervals")
      if (ntp > 1 && scatter <= 1) {
        trealigner.nt = ntp
      }
      if (scatter > 1) {
        trealigner.scatterCount = scatter
      }
      add(trealigner)

      indelrealigner.input_file :+= a
      if (know != null) {
        indelrealigner.known = know
      }
      indelrealigner.targetIntervals = trealigner.out
      indelrealigner.out = swapExt(outdir, a, "bam", "realn.bam")
      if (scatter > 1) {
        indelrealigner.scatterCount = scatter
      }
      add(indelrealigner)

      brecalibrator.input_file :+= indelrealigner.out
      if (sites != null) {
        brecalibrator.knownSites = sites
      }
      if (scatter > 1) {
        brecalibrator.scatterCount = scatter
      }
      brecalibrator.out = swapExt(outdir, a, "bam", "recal.cvs")
      add(brecalibrator)

      printreads.input_file :+= indelrealigner.out
      printreads.BQSR = brecalibrator.out
      printreads.out = swapExt(outdir, a, "bam", "recal.bam")
      if (scatter > 1) {
        printreads.scatterCount = scatter
      }
      add(printreads)

      if (!unifiedvcf) {
        val hapcaller = new HaplotypeCaller with referenceseq with hapcall //with nctparallel
        hapcaller.input_file :+= printreads.out
        hapcaller.out = swapExt(outdir, a, "bam", "HC.raw.vcf")
        if (scatter > 1) {
          hapcaller.scatterCount = scatter
        }
        add(hapcaller)
      }
      inputhap +:= printreads.out
    }

    if (unifiedvcf) {
      val hapcaller = new HaplotypeCaller with referenceseq with hapcall with intervals //with nctparallel
      hapcaller.input_file = inputhap
      hapcaller.out = outdir + "/" + hapcall
      if (scatter > 1) {
        hapcaller.scatterCount = scatter
      }
      add(hapcaller)
    }

    Runtime.getRuntime().exec("rm *.realn.bam")
    Runtime.getRuntime().exec("rm *.realn.bam.bai")
    Runtime.getRuntime().exec("rm *.recal.bam")
    Runtime.getRuntime().exec("rm *.recal.bam.bai")

  }
}
