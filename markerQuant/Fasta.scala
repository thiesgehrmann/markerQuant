// Written by Thies Gehrmann thiesgehrmann@gmail.com 2015-2016


import scala.util.parsing.combinator._
import scala.annotation._
import java.io._

package markerQuant{
object Fasta {


  case class Entry( description: String, sequence: BioSeq.DNASeq ) {
    def revcomp = new Entry(this.description, this.sequence.revcomp)
    override def equals(o: Any) = o match{
      case that: Fasta.Entry => this.sequence.equals(that.sequence)
      case _ => false
    }

    override def hashCode = this.sequence.hashCode

    override def toString = this.sequence.toString

  }

  /////////////////////////////////////////////////////////////////////////////

  class groupFastaLinesIterator(iter: BufferedIterator[String]) extends Iterator[Entry] {
    var nsequences = 0
    def hasNext = iter.hasNext
    def next = {
      if (!iter.hasNext) {
        Iterator.empty.next
      } else {

        @tailrec def untilNext(head: String, seq: String):  Entry = {
          val line = if(iter.hasNext) iter.head else ""
          if (!iter.hasNext || ( line.length > 1 && line.charAt(0) == '>')) {
            if (head.length == 0) {
              untilNext(iter.next.substring(1), "")
            } else {
              this.nsequences += 1
              if (nsequences % 5000 == 0){
                Utils.message("\rProcessed %d sequences".format(nsequences), ln=false)
              } else if (!iter.hasNext) {
                Utils.message("\rProcessed %d sequences".format(nsequences))
              }
              Entry(head, BioSeq.DNASeq.fromString(seq))
            }
          } else {
            untilNext(head, seq + iter.next)
          }
        }

        untilNext("", "")
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def read( fn: String ): groupFastaLinesIterator = {
    new groupFastaLinesIterator(io.Source.fromFile(fn).getLines.buffered)
  }

  /////////////////////////////////////////////////////////////////////////////

  def readMap(fn: String): Map[String,Entry] = {
    read(fn).map( e => (e.description, e)).toMap
  }

  /////////////////////////////////////////////////////////////////////////////

  def print(fastaList: Iterable[Entry]) = {
    def printSingleFasta(fastaSingle: Entry) = {
      println('>' + fastaSingle.description)
      println(fastaSingle.sequence.toString)
    }
    fastaList.map(printSingleFasta)
  }

  /////////////////////////////////////////////////////////////////////////////

  def write(fastaList: Iterable[Entry], fn: String) = {

    val outfd = new PrintWriter(new FileWriter(fn, false))

    fastaList.foreach{ e =>
      outfd.write(">" + e.description + '\n');
       e.sequence.toString.grouped(80).toList.foreach { l =>
        outfd.write(l + '\n');
      }
    }
    outfd.close()
  }

}
}
