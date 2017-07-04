import scala.annotation._

package markerQuant {

object Fastq {


  case class Entry( description: String, sequence: BioSeq.DNASeq, quality: Array[Byte] )

  /////////////////////////////////////////////////////////////////////////////

  class groupFastqLinesIterator(iter: BufferedIterator[String]) extends Iterator[Entry] {
    var nsequences = 0
    var nlines = 0
    def hasNext = iter.hasNext
    def next = {
      if (!iter.hasNext) {
        Iterator.empty.next
      } else {

        @tailrec def untilNext(seqid: String, seq: String, quality:String, state: Int): Fastq.Entry = {
          if (!iter.hasNext) {
            Iterator.empty.next
          } else {
            val line = if(iter.hasNext) iter.next else ""
            this.nlines += 1

            state match {
              case 0 => {
                if (line.charAt(0) == '@'){
                  untilNext(line, "", "", 1)
                } else {
                  untilNext("", "", "", 0)
                }
              }
              case 1 => {
                untilNext(seqid, line, "", 2)
              }
              case 2 => {
                if(line.charAt(0) == '+'){
                  untilNext(seqid, seq, "", 3)
                } else {
                  Utils.warning("Error on line %d. Expected '+' symbol.".format(this.nlines))
                  untilNext(seqid, seq, "", 2)
                }
              }
              case 3 => {
                this.nsequences += 1
                if (this.nsequences % 5000 == 0) { Utils.message("Processed %d fastq sequences".format(this.nsequences)) }
                Fastq.Entry(seqid, BioSeq.DNASeq.fromString(seq), quality.getBytes.map(x => (x - 33).toByte))
              }
            } 
          }
        }

        untilNext("", "", "", 0)
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def read( fn: String ): groupFastqLinesIterator = {
    Utils.message("Reading Fastq file '%s'".format(fn))
    new groupFastqLinesIterator(io.Source.fromFile(fn).getLines.buffered)
  }
  

}

}
