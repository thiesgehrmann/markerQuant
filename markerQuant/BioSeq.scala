package markerQuant {

object BioSeq {

    // I put this in place instead, to reduce the number of swaps between the spaces
  class DNASeq(val seq: String, val length: Int) {
    override def equals(o: Any) = o match{
      case that: DNASeq => this.seq.equals(that.seq)
      case _ => false
    }
   
    override def hashCode = this.seq.hashCode

    override def toString = this.seq

    def revcomp = new DNASeq(BioSeq.revcomp(this.seq), this.length)

  }

  object DNASeq {

    def apply(seq: String) = fromString(seq)

    def fromString(seq: String) = new DNASeq(seq.toLowerCase, seq.length)
  }

//  class DNASeq(val seq: Array[Int], val length: Int){
//
//    override def equals(o: Any) = o match {
//      case that: DNASeq => java.util.Arrays.equals(this.seq, that.seq) && this.length == that.length
//      case _ => false
//    }
//    override def hashCode = java.util.Arrays.hashCode(this.seq)
//
//    override def toString = {
//      this.seq.foldLeft(""){ (prev: String, next: Int) =>
//        //println(next.toBinaryString)
//        prev + (0 to DNASeq.nBasesPerVal-1).map{ i =>
//          val cbits = (next & (DNASeq.filledBits << (DNASeq.nBitsPerBase*i))) >>> (DNASeq.nBitsPerBase*i)
//          //println("%d: %s -> %s".format(i, cbits.toBinaryString, DNASeq.bit2base(cbits)))
//          DNASeq.bit2base(cbits)
//        }.mkString("")
//      }.slice(0,this.length)
//    }
//
//    /////////////////////////////////////////////////////////////////////////////
//
//    def revcomp = {
//      DNASeq.fromString(BioSeq.revcomp(this.toString))
//    }
//  }
//
//  /////////////////////////////////////////////////////////////////////////////
//
//  object DNASeq{
//
//    val nBitsPerBase = 3
//    val nBasesPerVal = 32/nBitsPerBase
//    val filledBits = Math.pow(2,nBitsPerBase).toInt-1
//
//    /////////////////////////////////////////////////////////////////////////////
//
//    def base2bit(base: Char) : Int = {
//      base match {
//        case 'a' => 0
//        case 'c' => 1
//        case 'g' => 2
//        case 't' => 3
//        case 'n' => 4
//        case '*' => 5
//        case _   => { Utils.error("The character '%c' is not a valid DNA base. Will be interpreted as 'n'.".format(base)); base2bit('n') }
//      }
//    }
//
//    ///////////////////////////////////////////////////////////////////////////
//
//    def bit2base(bit: Int) = {
//      bit match {
//        case 0 => 'a'
//        case 1 => 'c'
//        case 2 => 'g'
//        case 3 => 't'
//        case 4 => 'n'
//        case 5 => '*'
//        case _ => 'N'
//      }
//    }
//
//    ///////////////////////////////////////////////////////////////////////////
//
//    def fromString(seq: String) = {
//      val lseq = seq.toLowerCase
//      val bseq = (0 to seq.length by this.nBasesPerVal).map{ i =>
//        lseq.substring(i,math.min(i+this.nBasesPerVal, seq.length)).toCharArray.zipWithIndex.foldLeft(0){ (prev: Int, next: (Char,Int)) =>
//          //println("%c -> %d -> %32s | %32s".format(next._1, base2bit(next._1), prev.toBinaryString, (base2bit(next._1) << (next._2*2)).toBinaryString))
//          prev | (base2bit(next._1) << (next._2*this.nBitsPerBase))
//        }
//      }.toArray
//      new DNASeq(bseq, seq.length)
//    }
//
//  }

  ///////////////////////////////////////////////////////////////////////////

  def revcomp(seq: String) = {
    seq.reverse.toLowerCase.map{ b =>
      b match {
        case 'a' => 't'
        case 't' => 'a'
        case 'g' => 'c'
        case 'c' => 'g'
        case 'n' => 'n'
        case '*' => '*'
      }
    }
  }


}

}
