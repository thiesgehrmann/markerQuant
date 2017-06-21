package markerQuant {

object Test extends ActionObject {

  override def main(args: Array[String]) = {
    args(0) match {
      case "1" => { println("test1")
                    this.test1 }
      case opt => { println("Unknown test") }
    }
  }

  def test1 = {

    val sequence = Fasta.Entry("test", Fasta.FastaSeq.fromString("AAAAAAAAAATAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAT"))
    val kmers    = Markers.genKmers(sequence.sequence, 21)

    println(sequence.toString)
    kmers.foreach{ k => println("%d: %s".format(k.index, k.seq.toString)) }

    println("\nBut the set")
    kmers.toSet.foreach{ k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    println("\nAnd the Unique:")
    val uniqueKmers = Markers.uniqueKmers(kmers)
    uniqueKmers.foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    val kmerCounts = Markers.kmerCounts(Array(sequence, Fasta.Entry("test", Fasta.FastaSeq.fromString("TGATGTGATGGATAAAAAAAAAAAAAAAAAAAAATAGaaaaaaaaaaaaaaataaaaa"))), 21)
    println("\nAnd the counts:")
    kmerCounts.foreach{ case (k,v) => println("%d: %s: %d".format(k.index, k.seq.toString, v))}

    println("\nAnd the single count markers:")
    val single = Markers.singleCountKmers(uniqueKmers, kmerCounts)
    Markers.singleCountKmers(uniqueKmers, kmerCounts).foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    val aggregated = Markers.aggregateRedundant(single)
    println("\nAnd the aggregated markers")
    aggregated.foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    val gaps = aggregated.map(Markers.aggregatedKmerGaps(_, 21, kmerCounts)).flatten
    println("\nAnd the gaps to skip!")
    gaps.foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

  }

}

}
