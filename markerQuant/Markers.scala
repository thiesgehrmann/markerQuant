package markerQuant {

object Markers {

  case class Kmer(seq: Fasta.FastaSeq, index: Int) {
    /*************************************************************************
   *    * Define a kmertype such that we keep the index together with the kmer  *
   *       *   sequence. We have to redefine the equality and hashCode functions   *
   *          *   such that this type works with the sets and hashtables              *
   *             *************************************************************************/
    override def equals(o: Any) = o match {
      case that: Kmer => java.util.Arrays.equals(this.seq.seq, that.seq.seq) && this.seq.length == that.seq.length
      case _ => false
    }
    override def hashCode = java.util.Arrays.hashCode(this.seq.seq)
  }


  def genKmers(seq: Fasta.FastaSeq, k: Int) = {

    val upstream_downstream = (k-1)/2
    val startLocation = upstream_downstream
    val endLocation   = seq.length - upstream_downstream
    val stringSeq = seq.toString

    def genKmersHelper(index: Int): Stream[Kmer] = {
      if (index >= endLocation){
        Stream.empty
      } else {
        new Kmer(Fasta.FastaSeq.fromString(stringSeq.substring(index-upstream_downstream, index+upstream_downstream+1)), index) #:: genKmersHelper(index+1)
      }
    }

    genKmersHelper(startLocation)
  }

  def uniqueKmers(kmers: Seq[Kmer]) = {
    /* Returns an array of kmers that are unique in the array */

      // I feel like this can be optimized
    kmers.toSet.filter{item : Kmer => kmers.indexOf(item) == kmers.lastIndexOf(item)}.toArray.sortBy(_.index)
  }

  def kmerCounts(seqs: Seq[Fasta.Entry], k: Int): Map[Kmer,Int] = {
    seqs.map{ s => 
      genKmers(s.sequence, k).foldLeft(Map.empty[Kmer,Int]){ case (counts,kmer) =>
        counts ++ Map(kmer -> (counts.getOrElse(kmer,0)+1))
      }
    }.reduce( Utils.sumDicts(_, _))
  }

  def singleCountKmers(kmers: Seq[Kmer], kmerCounts: Map[Kmer,Int]) = {
    kmers.filter(kmerCounts.getOrElse(_,1) == 1)
  }

  def aggregateRedundant(kmers: Seq[Kmer]) = {
    val k = kmers(0).seq.length
    kmers.sortBy(_.index).foldLeft(List.empty[Kmer]){ case (l, kmer) =>
      if (l.length == 0) {
        List(kmer)
      } else {
        val prev = l.head
        val distance = kmer.index - prev.index
        if (distance < k) { // The kmers are overlapping
          val prevSeq = prev.seq.toString
          val kmerSeq = kmer.seq.toString
          val newSeq  = prevSeq + kmerSeq.substring(k-distance)
          new Kmer(Fasta.FastaSeq.fromString(newSeq), kmer.index) :: l.tail
        } else {
          kmer :: l
        }
      }
    }
  }

  def aggregatedKmerGaps(kmer: Kmer, k: Int, kmerCounts: Map[Kmer,Int]) = {
    val aggStartIndex = kmer.index - kmer.seq.length + ((k-1)/2)+1
    genKmers(kmer.seq, k).filter(kmerCounts.getOrElse(_,1) > 1).map( e => new Kmer(e.seq, e.index + aggStartIndex))
  }

}

}
