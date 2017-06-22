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
    }

    val fastaKmer = Fasta.read(options("markers")).toArray
    val fastaGaps = Fasta.read(options("gaps")).toArray
    val k         = options("k").toInt
    val sSpecific = if (options("strandSpecific") == "True") true else false
    val fastq     = options("fastq").split(',').map(Fastq.read)
    val outFile   = options("out")
    val counts    = new Utils.CountMap[String]

    val tree = Quantification.genTree(fastaKmer, fastaGaps, k, sSpecific)

    while(fastq.forall(_.hasNext)){
      val reads   = fastq.map(_.next).toArray
      val markers = Quantification.findInTree(tree, reads)
      counts.add(markers)
    }

    val outfd = new PrintWriter(new FileWriter(outFile, false))
    counts.get.foreach{ case (id, count) =>
      outfd.write("%s\t%d\n".format(id, count))
    }
    outfd.close()

  }

  val defaultArgs = Map("singleStrand" -> "False",
                        "out" -> "./quantification.tsv",
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
        case "-o" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("out" -> value))
        }
        case "-k" :: value :: tail => {
          processArgsHelper(tail, options ++ Map("k" -> value))
        }
        case "-s" :: tail => {
          processArgsHelper(tail, options ++ Map("singleStrand" -> "True"))
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

  }

}

}
