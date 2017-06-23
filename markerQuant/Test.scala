import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult

import scala.collection.JavaConverters._

package markerQuant {

object Test extends ActionObject {

  override def main(args: Array[String]) = {
    args(0) match {
      case "markergen" => this.markergen
      case "treetest"  => this.treetest
      case opt         => { println("Unknown test") }
    }
  }

  def markergen = {

    val sequence = Fasta.Entry("test", BioSeq.DNASeq.fromString("AAAAAAAAAATGAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAT"))
    val genome   = Array(sequence, Fasta.Entry("test2", BioSeq.DNASeq.fromString("TGATGTGATGGATAAAAAAAAAAAAAAAAAAAAATAGaaaaaaaaaaaaaaataaaaa")))
    val kmers    = Markers.genKmers(sequence.sequence, 21)

    println(sequence.toString)
    kmers.foreach{ k => println("%d: %s".format(k.index, k.seq.toString)) }

    println("\nBut the set")
    kmers.toSet.foreach{ k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    println("\nAnd the Unique:")
    val uniqueKmers = Markers.uniqueKmers(kmers).toArray.sortBy(_.index)
    uniqueKmers.foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    val kmerCounts = Markers.kmerCounts(genome ++ genome.map(_.revcomp), 21)
    println("\nAnd the counts:")
    kmerCounts.keys.toArray.sortBy(_.index).foreach{ k => println("%d: %s: %d".format(k.index, k.seq.toString, kmerCounts(k)))}

    println("\nAnd the single count markers:")
    val single = Markers.singleCountKmers(uniqueKmers, kmerCounts, false)
    single.foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    val aggregated = Markers.aggregateRedundant(single)
    println("\nAnd the aggregated markers")
    aggregated.foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

    val gaps = aggregated.map(Markers.aggregatedKmerGaps(_, 21, kmerCounts, false)).flatten
    println("\nAnd the gaps to skip!")
    gaps.foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

  }


  def treetest = {

    val tree = new AhoCorasick()
    tree.add("AAAAAA".getBytes, "id1")
    tree.add("TTTTTT".getBytes, "id1")

    tree.prepare()

    tree.search("AAAAAATTTTTT".getBytes).asInstanceOf[java.util.Iterator[SearchResult]].asScala.toArray.foreach( r => r.asInstanceOf[SearchResult].getOutputs.toArray.foreach{ rs => println(rs.asInstanceOf[String])})

  }

}

}
