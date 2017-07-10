// Written by Thies Gehrmann thiesgehrmann@gmail.com 2015-2016

/******************************************************************************
 * Determine Unique sequence identifiers based on on two inputs               *
 *                                                                            *
 ******************************************************************************/

import scala.annotation.tailrec

package markerQuant {

object UniqueMarkers extends ActionObject {

  override val description = "Find unique markers a set of genes and genomes"

  /////////////////////////////////////////////////////////////////////////////

  override def main(args: Array[String]) {

    val arguments = processArgs(args.toList)

    if (!arguments.contains("targets") || !arguments.contains("genomes")) {
      Utils.error("You need to define a targets file and a genome file")
      System.exit(1)
    }

      // Read input
    val inputFastas  = arguments("targets").split(',')
    val inputGenomes = arguments("genomes").split(',')
    val inputTrans   = arguments("transcriptomes").split(',')
    val k            = arguments("k").toInt
    val markerOut    = arguments("markerOut")
    val gapOut       = arguments("gapOut")
    val strSpecific  = if (arguments("strandSpecific") == "True") true else false

      // Read the target genes
    val targets  = inputFastas.map(Fasta.read(_).toArray).flatten
      // Read genome FASTA files
    val genomes = inputGenomes.map(Fasta.read(_).toArray).flatten;
      // Read the transcriptome FASTA files
    val transcriptomes = inputTrans.filter(_.length > 0).map(Fasta.read(_).toArray).flatten

    Utils.message("Counting Kmers.")
    val kmerCountsGenomes = Markers.kmerCounts(genomes ++ genomes.map(_.revcomp), k)
    val kmerCountsTrans   = if(strSpecific) {Markers.kmerCounts(transcriptomes, k)} else { Markers.kmerCounts(transcriptomes ++ transcriptomes.map(_.revcomp), k) }
    val kmerCountsTargets = if(strSpecific) {Markers.kmerCounts(targets, k)} else {Markers.kmerCounts(targets ++ targets.map(_.revcomp), k) }

    val kmerCounts = Array(kmerCountsTargets, kmerCountsGenomes, kmerCountsTrans).filter(_.size > 0)

    val (targetKmers, gapKmers) = Markers.getUniqueAggregatedKmersAndGaps(targets, kmerCounts, k, strSpecific)

    //Utils.message("Filtering generated Kmers for the target sequences.")
    //val filteredKmers = targets.map( s => kmerCounts.foldLeft(Markers.genKmers(s.sequence, k): Iterable[Markers.Kmer]){ case (kmers, counts) => Markers.singleCountKmers(kmers, counts, strSpecific)})
    //val filteredKmers = targets.map(s => Markers.singleCountKmers(Markers.singleCountKmers(Markers.uniqueKmers(Markers.genKmers(s.sequence, k)), kmerCountsGenomes, strSpecific), kmerCountsTargets, strSpecific))
    //filteredKmers.foreach(group => group.foreach{ k => println("%d: %s".format(k.index, k.seq.toString))})

    //Utils.message("Aggregating the Kmers.")
    //val targetKmers = filteredKmers.map(Markers.aggregateRedundant(_))

    //Utils.message("Finding gaps in the aggregated kmers")
    //val gapKmers = targetKmers.map( group => group.map( kmer => kmerCounts.map( counts => Markers.aggregatedKmerGaps(kmer, k, counts, strSpecific).toSet).reduce(_ union _)))
    //val gapKmers = targetKmers.map( group => group.map( kmer => Markers.aggregatedKmerGaps(kmer, k, kmerCountsGenomes, strSpecific).toSet union Markers.aggregatedKmerGaps(kmer, k, kmerCountsTargets, strSpecific).toSet))

    //println("Kmers to use:")
    val faKmers = targetKmers.zipWithIndex.map{ case (kmers, index) =>
      val targetDescription = targets(index).description
      kmers.map{ kmer => Fasta.Entry("%s:%d".format(targetDescription, kmer.index), kmer.seq) }
    }.flatten
    Fasta.write(faKmers, markerOut)

    //targetKmers.flatten.map{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    //println("Gap kmers to skip:")
    val faGaps = gapKmers.zipWithIndex.map{ case (kmers, index) =>
      kmers.map{ gap =>
        Fasta.Entry("%s:%d".format(targets(index).description, gap.index), gap.seq)
      }
    }.flatten
    Fasta.write(faGaps, gapOut)
    //faGaps.foreach( g => println("%s: %s".format(g.description, g.sequence)))

    System.exit(0);

  }

  /////////////////////////////////////////////////////////////////////////////

  val defaultArgs = Map("transcriptomes" -> "",
                        "strandSpecific" -> "False",
                        "markerOut" -> "./markers.fasta",
                        "gapOut" -> "./gaps.fasta",
                        "k" -> "21")

  def processArgs(args: List[String]) : Map[String,String] = {

    @tailrec def processArgsHelper(args: List[String], options: Map[String,String]): Map[String,String] = {

      args match {
        case "-t" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("targets" -> value))
        }
        case "-T" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("transcriptomes" -> value))
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
          processArgsHelper(tail, options ++ Map("strandSpecific" -> "True"))
        }
        case "-M" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("markerOut" -> value))
        }
        case "-G" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("gapOut" -> value))
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
    println("  -T <fastaFile>: A fasta file containing the transcriptome sequences you want the markers to be unique for. [None]")
    println("  -k <fastaFile>: Size of kmer to use (must be odd) [%s]".format(defaultArgs("k")))
    println("  -M <outFile>: Where to output the markers [%s]".format(defaultArgs("markerOut")))
    println("  -G <outFile>: Where to output the gaps [%s]".format(defaultArgs("gapOut")))
    println("  -s :            Generate strand specific markers.")
    println("")
    println("Example usage: quant -t mytargets.fasta -g mygenome.fasta # Generates unique markers relative to the genome")
    println("               quant -t mytargets.fasta -g mygenome.fasta -T mytranscriptome.fasta # Generates unique markers relative to genome,transcriptome")
    println("               quant -t mytargets.fasta -g mygenome.fasta -M m.fasta -G g.fasta -s # Generates unique markers for strand specific RNA-Seq, but puts them in different files.")

  }

  /////////////////////////////////////////////////////////////////////////////*/

}

}
