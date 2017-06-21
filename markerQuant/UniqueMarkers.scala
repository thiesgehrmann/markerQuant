// Written by Thies Gehrmann thiesgehrmann@gmail.com 2015-2016

/******************************************************************************
 * Determine Unique sequence identifiers based on on two inputs               *
 *                                                                            *
 ******************************************************************************/

import scala.annotation.tailrec

package markerQuant {

object UniqueMarkers extends ActionObject {

  override def description = "Find unique markers a set of genes and genomes"

  /////////////////////////////////////////////////////////////////////////////

  override def main(args: Array[String]) {

    val arguments = processArgs(args.toList)

    if (!arguments.contains("targets") || !arguments.contains("genomes")) {
      Utils.error("You need to define a targets file and a genome file")
    }

      // Read input
    val inputFastas  = arguments("targets").split(',')
    val inputGenomes = arguments("genomes").split(',')
    val k            = arguments("k").toInt
    val outPrefix    = arguments("outprefix")
    val strSpecific  = if (arguments("singleStrand") == "True") true else false

      // Read the target genes
    val targets  = inputFastas.map(Fasta.read(_).toArray).flatten
      // Read genome FASTA files
    val genomes = inputGenomes.map(Fasta.read(_).toArray).flatten;

    val kmerCounts  = Markers.kmerCounts(genomes, k)
    val targetKmers = targets.map(s => Markers.aggregateRedundant(Markers.singleCountKmers(Markers.uniqueKmers(Markers.genKmers(s.sequence, k)), kmerCounts)))
    val gapKmers    = targetKmers.flatten.map( Markers.aggregatedKmerGaps(_, k, kmerCounts)).flatten

    println("Kmers to use:")
    val faKmers = targetKmers.zipWithIndex.map{ case (kmers, index) =>
      val targetDescription = targets(index).description
      kmers.map{ kmer => Fasta.Entry("%s:%d".format(targetDescription, kmer.index), kmer.seq) }
    }.flatten
    Fasta.write(faKmers, "%s%s".format(outPrefix, "markers.fasta"))

    targetKmers.flatten.map{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    println("Gap kmers to skip:")
    val faGaps = gapKmers.zipWithIndex.map{ case (kmer, index) =>
      val targetDescription = targets(index).description
      Fasta.Entry("%s:%d".format(targetDescription, kmer.index), kmer.seq)
    }
    Fasta.write(faGaps, "%s%s".format(outPrefix, "gaps.fasta"))
    gapKmers.map{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    System.exit(0);

  }

  /////////////////////////////////////////////////////////////////////////////

  val defaultArgs = Map("singleStrand" -> "False",
                        "outprefix" -> "./",
                        "k" -> "21")

  def processArgs(args: List[String]) : Map[String,String] = {

    @tailrec def processArgsHelper(args: List[String], options: Map[String,String]): Map[String,String] = {

      args match {
        case "-t" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("targets" -> value))
        }
        case "-g" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("genomes" -> value))
        }
        case "-k" :: value ::tail => {
          val size = value.toInt
          if (size % 2 == 0 || size < 3) {
            Utils.error("'-k' must have odd value and be larger than 3!")
          }
          processArgsHelper(tail, options ++ Map("k" -> value))
        }
        case "-s" :: tail => {
          processArgsHelper(tail, options ++ Map("singleStrand" -> "True"))
        }
        case "-o" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("outprefix" -> value))
        }
        case opt :: tail => {
          Utils.error("Unknown option '%s'.".format(opt))
          processArgsHelper(tail, options)
        }
        case Nil => { options  }
      }
    }

    processArgsHelper(args, this.defaultArgs)

  }

  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////*/

  override def usage = {
    println("Unique identifier detection")
    println("Usage: markers <options>");
    println("  -t <fastaFile>: A fasta file containing the gene sequences you want to have markers for (multiple can be provided, separated by commas).")
    println("  -g <fastaFile>: A fasta file containing the genome sequences you want the markers to be unique for.")
    println("  -k <fastaFile>: Size of kmer to use [21]")
    println("  -o <fastaFile>: Where to output the files [./]")
    println("  -s :            Generate strand specific markers.")
    println("")

  }

  /////////////////////////////////////////////////////////////////////////////*/

}

}
