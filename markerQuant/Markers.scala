package markerQuant {

object Markers {

  case class Kmer(seq: BioSeq.DNASeq, index: Int) {
    /*************************************************************************
   *    * Define a kmertype such that we keep the index together with the kmer  *
   *       *   sequence. We have to redefine the equality and hashCode functions   *
   *          *   such that this type works with the sets and hashtables              *
   *             *************************************************************************/
    override def equals(o: Any) = {
      o match{
        case that: Kmer => this.seq.equals(that.seq)
        case _ => false
      }
    }
    override def hashCode = this.seq.hashCode

    def revcomp = {
      new Kmer(this.seq.revcomp, this.index)
    }

    def toFastaEntry = {
      Fasta.Entry("%d".format(this.index), this.seq)
    }
    def toFastqEntry = {
      this.toFastaEntry.toFastqEntry
    }
  }

  def genKmers(seq: BioSeq.DNASeq, k: Int) = {

    val upstream_downstream = (k-1)/2
    val startLocation = upstream_downstream
    val endLocation   = seq.length - upstream_downstream
    val stringSeq = seq.toString

    def genKmersHelper(index: Int): Stream[Kmer] = {
      if (index >= endLocation){
        Stream.empty
      } else {
        // Skip 'n' characters
        val substr = stringSeq.substring(index-upstream_downstream, index+upstream_downstream+1)
        val nIndex = substr.indexOf('n')
        if (nIndex > 0) {
          genKmersHelper(index+nIndex+upstream_downstream)
        } else {
          new Kmer(BioSeq.DNASeq.fromString(substr), index) #:: genKmersHelper(index+1)
        }
      }
    }

    genKmersHelper(startLocation)
  }

  def uniqueKmerGroups(kmers: Seq[Seq[Kmer]]) : Seq[Seq[Kmer]] = {
    val counts = new Utils.CountMap[Kmer]
    kmers.foreach(group => counts.add(group))
    kmers.map(group => group.filter(kmer => (counts(kmer) <= 1)))
  }

  def uniqueKmers(kmers: Seq[Kmer]) = {
    /* Returns an array of kmers that are unique in the array */

      // I feel like this can be optimized
    kmers.toSet.filter{item : Kmer => kmers.indexOf(item) == kmers.lastIndexOf(item)}.toArray.sortBy(_.index)
  }

  def kmerCounts(seqs: Iterable[Fasta.Entry], k: Int): Map[Kmer,Int] = {
    if (seqs.isEmpty) {
      Map.empty[Kmer,Int]
    } else {
      seqs.zipWithIndex.map{ case (s,i) => 
        genKmers(s.sequence, k).foldLeft(Map.empty[Kmer,Int]){ case (counts,kmer) =>
          counts ++ Map(kmer -> (counts.getOrElse(kmer,0)+1))
        }
      }.reduce( Utils.sumDicts(_, _))
    }
  }

  

  def singleCountKmers(kmers: Iterable[Kmer], kmerCounts: Map[Kmer,Int], strandSpecific: Boolean) = {
    if (strandSpecific) {
      kmers.filter{ kmer =>
        //println("%d:%s -> %d".format(kmer.index, kmer.seq.toString, kmerCounts.getOrElse(kmer,1)))
        kmerCounts.getOrElse(kmer,0) <= 1}
    } else {
      kmers.filter{kmer =>
        //println("%d: %s -> %d".format(kmer.index, kmer.seq.toString, kmerCounts.getOrElse(kmer,0) + kmerCounts.getOrElse(kmer.revcomp,0)))
        (kmerCounts.getOrElse(kmer,0) <= 1) && (kmerCounts.getOrElse(kmer.revcomp,0) <= 1) }
    }
  }

  def aggregateRedundant(kmers: Iterable[Kmer]) = {
    val kmerArray = kmers.toArray
    val k = if (kmerArray.length == 0) 21 else kmerArray(0).seq.length
    kmerArray.sortBy(_.index).foldLeft(List.empty[Kmer]){ case (l, kmer) =>
      if (l.length == 0) {
        List(kmer)
      } else {
        val prev = l.head
        val distance = kmer.index - prev.index
        if (distance < k) { // The kmers are overlapping
          val prevSeq = prev.seq.toString
          val kmerSeq = kmer.seq.toString
          val newSeq  = prevSeq + kmerSeq.substring(k-distance)
          new Kmer(BioSeq.DNASeq.fromString(newSeq), kmer.index) :: l.tail
        } else {
          kmer :: l
        }
      }
    }
  }

  def aggregatedKmerGaps(kmer: Kmer, k: Int, kmerCounts: Array[Map[Kmer,Int]]) = {
    val aggStartIndex = kmer.index - kmer.seq.length + ((k-1)/2)+1
    genKmers(kmer.seq, k).filter{ kmer => kmerCounts.map( kc => kc.getOrElse(kmer,1) > 1).foldLeft(false)(_ || _)}.map( e => new Kmer(e.seq, e.index + aggStartIndex)).toSet
  }

  def getUniqueAggregatedKmersAndGaps(targets: Iterable[Fasta.Entry], kmerCountsArg: Array[Map[Markers.Kmer,Int]], k: Int, strSpecific: Boolean) = {

    val kmerCounts = kmerCountsArg.filter(_.size > 0)

    val filteredKmers = targets.map( s => kmerCounts.foldLeft(Markers.genKmers(s.sequence, k): Iterable[Markers.Kmer]){ case (kmers, counts) => Markers.singleCountKmers(kmers, counts, strSpecific)})
    
    Utils.message("Aggregating the Kmers.")
    val targetKmers = filteredKmers.map(Markers.aggregateRedundant(_))

    Utils.message("Finding gaps in the aggregated kmers")
    val gapKmers = targetKmers.map( group => group.map( kmer => aggregatedKmerGaps(kmer, k, kmerCounts).toSet).foldLeft(Set.empty[Kmer])(_ union _))

    (targetKmers, gapKmers)

  }

}

}
