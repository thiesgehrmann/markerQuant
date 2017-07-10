import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult

import scala.annotation._
import java.io._

package markerQuant {

object Quantify extends ActionObject {

  override val description = "Quantify markers in a fastq file"

  override def main(args: Array[String]) = {
    val options = processArgs(args.toList)

    if (!(( options contains "markers") && (options contains "gaps") && (options contains "fastq"))){
      Utils.error("You must define -m, -g and -f")
      System.exit(1)
    }

    val fastaKmer     = Fasta.read(options("markers")).toArray
    val fastaGaps     = Fasta.read(options("gaps")).toArray
    val k             = options("k").toInt
    val sSpecific     = if (options("strandSpecific") == "True") true else false
    val fastq         = options("fastq").split(',').map(Fastq.read)
    val markerOutFile = options("markerOut")
    val targetOutFile = options("targetOut")
    val minQual       = options("minQual").toInt.toByte
    val mCounts       = Utils.CountMap[String]
    val tCounts       = Utils.CountMap[String]

    val quantifier = Quantification.Quantifier(fastaKmer, fastaGaps, k, sSpecific, fastq.length > 1, minQual)

    Utils.message("Counting markers")
    while(fastq.forall(_.hasNext)){
      val reads   = fastq.map(_.next).toArray
      val markers = quantifier.search(reads)
      val targets = markers.map(_.split(':')(0)).distinct

      mCounts.add(markers)
      //println("markers: %s".format(markers.mkString(", ")))
      tCounts.add(targets)
      //println("targets: %s -> %d".format(targets.mkString(", "), tCounts(if (targets.size > 0) { targets(0)} else "")))
    }

    /* Write the marker counts out */
    val mOutfd = new PrintWriter(new FileWriter(markerOutFile, false))
    fastaKmer.map(_.description).foreach{ id =>
      mOutfd.write("%s\t%d\n".format(id, mCounts(id)))
    }
    mOutfd.write("%s\t%d\n".format(Quantification.Quantifier.multiTargetRead, mCounts(Quantification.Quantifier.multiTargetRead)))
    mOutfd.write("%s\t%d\n".format(Quantification.Quantifier.unmappedRead, mCounts(Quantification.Quantifier.unmappedRead)))
    mOutfd.write("%s\t%d\n".format(Quantification.Quantifier.lowQualRead, mCounts(Quantification.Quantifier.lowQualRead)))
    mOutfd.close()

    /* Write the target counts out */
    val tOutfd = new PrintWriter(new FileWriter(targetOutFile, false))
    fastaKmer.map(f => f.description.substring(0, f.description.indexOf(':'))).distinct.foreach{ id =>
      tOutfd.write("%s\t%d\n".format(id, tCounts(id)))
    }
    tOutfd.write("%s\t%d\n".format(Quantification.Quantifier.multiTargetReadTarget, tCounts(Quantification.Quantifier.multiTargetReadTarget)))
    tOutfd.write("%s\t%d\n".format(Quantification.Quantifier.unmappedReadTarget, tCounts(Quantification.Quantifier.unmappedReadTarget)))
    tOutfd.write("%s\t%d\n".format(Quantification.Quantifier.lowQualReadTarget, tCounts(Quantification.Quantifier.lowQualReadTarget)))
    tOutfd.close()

  }

  val defaultArgs = Map("strandSpecific" -> "False",
                        "targetOut" -> "./target_quantification.tsv",
                        "markerOut" -> "./marker_quantification.tsv",
                        "minQual" -> "25",
                        "k" -> "21")

  def processArgs(args: List[String]) : Map[String,String] = {

    @tailrec def processArgsHelper(args: List[String], options: Map[String,String]): Map[String,String] = {

      args match {
        case "-m" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("markers" -> value))
        }
        case "-g" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("gaps" -> value))
        }
        case "-f" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("fastq" -> value))
        }
        case "-M" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("markerOut" -> value))
        }
        case "-T" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("targetOut" -> value))
        }
        case "-k" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("k" -> value))
        }
        case "-s" :: tail => {
          processArgsHelper(tail, options ++ Map("strandSpecific" -> "True"))
        }
        case "-q" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("minQual" -> value))
        }
        case opt :: tail => {
          Utils.warning("Option '%s' unknown.".format(opt))
          processArgsHelper(tail, options)
        }
        case Nil => { options }
      }

    }

    processArgsHelper(args, defaultArgs)

  }  

  override def usage = {
    println("Quantify unique markers")
    println("Usage: quant  <options>");
    println("  -m <markerFile>: A fasta file containing the gene sequences you want to have markers for (multiple can be provided, separated by commas).")
    println("  -g <gapFile>: A fasta file containing the genome sequences you want the markers to be unique for.")
    println("  -f <fastqFile1>[,fastqFile2]: A fasta file containing the transcriptome sequences you want the markers to be unique for. [None]")
    println("  -M <outFile>: Where to output counts for markers. [Default %s]".format(defaultArgs("markerOut")))
    println("  -T <outFile>: Where to output counts for targets. [Default %s]".format(defaultArgs("targetOut")))
    println("  -k <integer>: kmer size to use [Default %s]".format(defaultArgs("k")))
    println("  -s : Provided markers are strand-specific.")
    println("  -q <integer [0,42]>: Minimum quality the marker-mapped read region must have. [Default %s]".format(defaultArgs("minQual")))
    println("")

  }

}

}
