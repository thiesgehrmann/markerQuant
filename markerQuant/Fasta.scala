// Written by Thies Gehrmann thiesgehrmann@gmail.com 2015-2016


import scala.util.parsing.combinator._
import scala.annotation._
import java.io._

package markerQuant{
object Fasta {

  class FastaSeq(val seq: Array[Int], val length: Int){

    override def toString = {
      this.seq.foldLeft(""){ (prev: String, next: Int) =>
        //println(next.toBinaryString)
        prev + (0 to FastaSeq.nBasesPerVal-1).map{ i =>
          val cbits = (next & (FastaSeq.filledBits << (FastaSeq.nBitsPerBase*i))) >>> (FastaSeq.nBitsPerBase*i)
          //println("%d: %s -> %s".format(i, cbits.toBinaryString, FastaSeq.bit2base(cbits)))
          FastaSeq.bit2base(cbits)
        }.mkString("")
      }.slice(0,this.length)
    }

    /////////////////////////////////////////////////////////////////////////////

    def revcomp = {
      FastaSeq.fromString(Fasta.revcomp(this.toString))
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  object FastaSeq{

    val nBitsPerBase = 3
    val nBasesPerVal = 32/nBitsPerBase
    val filledBits = Math.pow(2,nBitsPerBase).toInt-1

    /////////////////////////////////////////////////////////////////////////////

    def base2bit(base: Char) : Int = {
      base match {
        case 'a' => 0
        case 'c' => 1
        case 'g' => 2
        case 't' => 3
        case 'n' => 4
        case '*' => 5
        case _   => { Utils.error("The character '%c' is not a valid DNA base. Will be interpreted as 'n'.".format(base)); base2bit('n') }
      }
    }

    ///////////////////////////////////////////////////////////////////////////

    def bit2base(bit: Int) = {
      bit match {
        case 0 => 'a'
        case 1 => 'c'
        case 2 => 'g'
        case 3 => 't'
        case 4 => 'n'
        case 5 => '*'
        case _ => 'N'
      }
    }

    ///////////////////////////////////////////////////////////////////////////

    def fromString(seq: String) = {
      val lseq = seq.toLowerCase
      val bseq = (0 to seq.length by this.nBasesPerVal).map{ i =>
        lseq.substring(i,math.min(i+this.nBasesPerVal, seq.length)).toCharArray.zipWithIndex.foldLeft(0){ (prev: Int, next: (Char,Int)) =>
          //println("%c -> %d -> %32s | %32s".format(next._1, base2bit(next._1), prev.toBinaryString, (base2bit(next._1) << (next._2*2)).toBinaryString))
          prev | (base2bit(next._1) << (next._2*this.nBitsPerBase))
        }
      }.toArray
      new FastaSeq(bseq, seq.length)
    }

  }

  /////////////////////////////////////////////////////////////////////////////

  case class Entry( description: String, sequence: FastaSeq )

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
              Entry(head, FastaSeq.fromString(seq))
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

  def revcomp(seq: String) = {
    seq.reverse.toLowerCase.map{ b =>
      b match {
        case 'a' => 'c'
        case 'c' => 'a'
        case 'g' => 't'
        case 't' => 'g'
      }
    }
  }

}
}
