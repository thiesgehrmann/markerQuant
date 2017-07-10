import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult

import scala.collection.JavaConverters._

package markerQuant {

object Test extends ActionObject {

  override def main(args: Array[String]) = {
    args(0) match {
      case "markergen" => this.markergen
      case "treetest"  => this.treetest
      case "markerset" => this.markerset
      case "markerquant" => this.markerquantTest
      case "uniquemarkertest" => this.uniqueMarkerTest
      case opt         => { println("Unknown test") }
    }
  }

  val minQual = 21.toByte

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

    val gaps = aggregated.map(Markers.aggregatedKmerGaps(_, 21, Array(kmerCounts))).flatten
    println("\nAnd the gaps to skip!")
    gaps.foreach{k: Markers.Kmer => println("%d: %s".format(k.index, k.seq.toString)) }

  }

  def markerset = {

    val reads = Array(Fastq.Entry("R1", BioSeq.DNASeq.fromString("gccccaatgttcgtcatgggtgttaacgaagaaaaatacacttctgacttgaagattgtctccaacgcttcttgtaccaccaaatgtttggctccattgg"), "".getBytes),
                      Fastq.Entry("R2", BioSeq.DNASeq.fromString("cacctctccagtccttgtgggatggaccgtcaacagtcttttgggtggcggtcatggagtgaacagtggtcatcaaaccttcttcaataccgaaagcatc"), "".getBytes))
    val markers = Array(Fasta.Entry("TDH1:444", BioSeq.DNASeq.fromString("gtcatcactgctccatcttcttctgctccaatgtttgttgttggtgttaaccacactaaatacactccagacaagaagattgtctccaacgcttcttgtaccac")),
                        Fasta.Entry("TDH2:555", BioSeq.DNASeq.fromString("gttcactccatgaccgccacccaaaagactgttgacggtcc")))

    val kmerCounts = Markers.kmerCounts(markers, 21)
    val (targetKmers, gapKmers) =  Markers.getUniqueAggregatedKmersAndGaps(markers, Array(kmerCounts), 21, false)

    targetKmers.flatten.foreach( k => println(">%d\n%s".format(k.index, k.seq.toString)))

    val targetKmerFasta = targetKmers.zipWithIndex.map{ case (group, i) => group.map( kmer => Fasta.Entry("%s:%d".format(markers(i).description, kmer.index), kmer.seq))}.flatten
    val gapKmerFasta = gapKmers.flatten.toSet.toArray.map(_.toFastaEntry)
    val quantifier = Quantification.Quantifier(targetKmerFasta, gapKmerFasta, 21, false, reads.length > 1, minQual, Array.empty[String].toSet)

    reads.foreach( read => quantifier.search(Array(read)))

    //reads.map(x => Markers.genKmers(x.sequence, 21)).flatten.toSet diff gapKmers.toSet

  }

  def markerquantTest = {

    val genes = Array(Fasta.Entry(">YJL052W_mRNA gene=TDH1", BioSeq.DNASeq("ATGATCAGAATTGCTATTAACGGTTTCGGTAGAATCGGTAGATTGGTCTTGAGATTGGCTTTGCAAAGAAAAGACATTGAGGTTGTTGCTGTCAACGATCCATTTATCTCTAACGATTATGCTGCTTACATGGTCAAGTACGATTCTACTCATGGTAGATACAAGGGTACTGTTTCCCATGACGACAAGCACATCATCATTGATGGTGTCAAGATCGCTACCTACCAAGAAAGAGACCCAGCTAACTTGCCATGGGGTTCTCTAAAGATCGATGTCGCTGTTGACTCCACTGGTGTTTTCAAGGAATTGGACACCGCTCAAAAGCACATTGACGCTGGTGCCAAGAAGGTTGTCATCACTGCTCCATCTTCTTCTGCTCCAATGTTTGTTGTTGGTGTTAACCACACTAAATACACTCCAGACAAGAAGATTGTCTCCAACGCTTCTTGTACCACCAACTGTTTGGCTCCATTGGCCAAGGTTATCAACGATGCTTTCGGTATTGAAGAAGGTTTGATGACCACTGTTCACTCCATGACCGCCACTCAAAAGACTGTTGATGGTCCATCCCACAAGGACTGGAGAGGTGGTAGAACCGCTTCCGGTAACATTATCCCATCCTCTACCGGTGCTGCTAAGGCTGTCGGTAAGGTCTTGCCAGAATTGCAAGGTAAGTTGACCGGTATGGCTTTCAGAGTCCCAACCGTCGATGTTTCCGTTGTTGACTTGACTGTCAAGTTGGAAAAGGAAGCTACTTACGACCAAATCAAGAAGGCTGTTAAGGCTGCCGCTGAAGGTCCAATGAAGGGTGTTTTGGGTTACACCGAAGATGCCGTTGTCTCCTCTGATTTCTTGGGTGACACTCACGCTTCCATCTTCGATGCCTCCGCTGGTATCCAATTGTCTCCAAAGTTCGTCAAGTTGATTTCCTGGTACGATAACGAATACGGTTACTCCGCCAGAGTTGTTGACTTGATCGAATATGTTGCCAAGGCTTAA")),
                      Fasta.Entry(">YJR009C_mRNA gene=TDH2", BioSeq.DNASeq("ATGGTTAGAGTTGCTATTAACGGTTTCGGTAGAATCGGTAGATTGGTTATGAGAATTGCTTTGCAAAGAAAGAACGTCGAAGTTGTTGCTTTGAACGATCCTTTCATCTCTAACGACTACTCCGCTTACATGTTCAAGTACGACTCTACTCACGGTAGATACGCTGGTGAAGTTTCCCACGATGACAAGCACATCATCGTTGATGGTCACAAGATCGCCACTTTCCAAGAAAGAGACCCAGCTAACTTGCCATGGGCTTCTCTAAACATTGACATCGCCATTGACTCCACTGGTGTTTTCAAGGAATTGGACACTGCTCAAAAGCACATTGACGCTGGTGCCAAGAAGGTTGTCATCACTGCTCCATCTTCCACCGCCCCAATGTTCGTCATGGGTGTTAACGAAGAAAAATACACTTCTGACTTGAAGATTGTTTCCAACGCTTCTTGTACCACCAACTGTTTGGCTCCATTGGCCAAGGTTATCAACGATGCTTTCGGTATTGAAGAAGGTTTGATGACCACTGTTCACTCCATGACCGCCACCCAAAAGACTGTTGACGGTCCATCCCACAAGGACTGGAGAGGTGGTAGAACCGCTTCCGGTAACATCATCCCATCCTCTACCGGTGCTGCTAAGGCTGTCGGTAAGGTCTTGCCAGAATTGCAAGGTAAGTTGACCGGTATGGCTTTCAGAGTCCCAACCGTCGATGTTTCCGTTGTTGACTTGACTGTCAAGTTGAACAAGGAAACCACCTACGATGAAATCAAGAAGGTTGTCAAGGCTGCCGCTGAAGGTAAGTTGAAGGGTGTCTTGGGTTACACTGAAGACGCTGTTGTCTCCTCTGACTTCTTGGGTGACTCTAACTCTTCCATCTTCGATGCTGCCGCTGGTATCCAATTGTCTCCAAAGTTCGTCAAGTTGGTTTCCTGGTACGACAACGAATACGGTTACTCTACCAGAGTTGTCGACTTGGTTGAACACGTTGCCAAGGCTTAA")))

    val kmerCounts = Markers.kmerCounts(genes, 21)
    val (markers, gaps) = Markers.getUniqueAggregatedKmersAndGaps(genes, Array(kmerCounts), 21, false)
    
    val quantifier = Quantification.Quantifier(markers.zipWithIndex.map{case (g, i) => g.map( k => Fasta.Entry("%s:%d".format(genes(i).description, k.index), k.seq))}.flatten, gaps.toArray.flatten.map(_.toFastaEntry), 21, false, false, minQual, Array.empty[String].toSet)

    genes.foreach{ g =>
      println("Checking kmers for %s".format(g.description))
      Markers.genKmers(g.sequence, 21).foreach{ k =>
        println("F-%4d: %s -> %s".format(k.index, k.seq.toString, quantifier.search(Array(k.toFastqEntry)).mkString(",")))
        println("R-%4d: %s -> %s".format(k.index, k.seq.toString, quantifier.search(Array(k.revcomp.toFastqEntry)).mkString(",")))
      }
    }
  }

  def treetest = {

    val tree = new AhoCorasick()
    tree.add("AAAAAA".getBytes, "id1")
    tree.add("TTTTTT".getBytes, "id1")

    tree.prepare()

    tree.search("AAAAAATTTTTT".getBytes).asInstanceOf[java.util.Iterator[SearchResult]].asScala.toArray.foreach( r => r.asInstanceOf[SearchResult].getOutputs.toArray.foreach{ rs => println(rs.asInstanceOf[String])})

  }

  def uniqueMarkerTest = {
    val genes = Array(Fasta.Entry("ADH1", BioSeq.DNASeq("ATGTCTATCCCAGAAACTCAAAAAGGTGTTATCTTCTACGAATCCCACGGTAAGTTGGAATACAAAGATATTCCAGTTCCAAAGCCAAAGGCCAACGAATTGTTGATCAACGTTAAATACTCTGGTGTCTGTCACACTGACTTGCACGCTTGGCACGGTGACTGGCCATTGCCAGTTAAGCTACCATTAGTCGGTGGTCACGAAGGTGCCGGTGTCGTTGTCGGCATGGGTGAAAACGTTAAGGGCTGGAAGATCGGTGACTACGCCGGTATCAAATGGTTGAACGGTTCTTGTATGGCCTGTGAATACTGTGAATTGGGTAACGAATCCAACTGTCCTCACGCTGACTTGTCTGGTTACACCCATGACGGTTCTTTCCAACAATACGCTACTGCTGACGCGGTGCAAGCCGCTCGTATTCCCGAAGGGACCGACTTGGCCCAAGTCGCCCCCATCTTGTGTGCTGGTATCACCGTCTACAAGGCTTTGAAGTCTGCTAACTTGATGGCCGGTCACTGGGTTGCTATCTCCGGTGCTGCTGGTGGTCTAGGTTCTTTGGCTGTTCAATACGCCAAGGCTATGGGTTACAGAGTCTTGGGTATTGACGGTGGTGAAGGTAAGGAAGAATTATTCAGATCCATCGGTGGTGAAGTCTTCATTGACTTCACTAAGGAAAAGGACATTGTCGGTGCTGTTCTAAAGGCCACTGACGGTGGTGCTCACGGTGTCATCAACGTTTCCGTTTCCGAAGCCGCTATTGAAGCTTCTACCAGATACGTTAGAGCTAACGGTACCACCGTTTTGGTCGGTATGCCAGCTGGTGCCAAGTGTTGTTCTGATGTCTTCAACCAAGTCGTCAAGTCCATCTCTATTGTTGGTTCTTACGTCGGTAACAGAGCTGACACCAGAGAAGCTTTGGACTTCTTCGCCAGAGGTTTGGTCAAGTCTCCAATCAAGGTTGTCGGCTTGTCTACCTTGCCAGAAATTTACGAAAAGATGGAAAAGGGTCAAATCGTTGGTAGATACGTTGTTGACACTTCTAAATAA")),
                      Fasta.Entry("ADH1_blocks", BioSeq.DNASeq("ATGTCTATCCCAGAAACTCAAAAAGGTGTTATCTTCTACGAATCCCACGGTAAGTTGGAATACAAAGATATTCCAGTTCCAAAGCCAAAGGCCAACGAATTGTTAATCAACGTTAAATACTCTGGTGTCTGTCACACTGACTTGCACGCTTGGCACGGTGACTGGCCATTACCAGTTAAGCTACCATTAGTCGGTGGTCACGAAGGTGCAGGTGTCGTTGTCGGCATGGGTGAAAACGTTAAGGGCTGGAAGATCGGTGACTACGCAGGTATCAAATGGTTGAACGGTTCTTGTATGGCCTGTGAATACTGTGAATTAGGTAACGAATCCAACTGTCCTCACGCTGACTTGTCTGGTTACACACATGACGGTTCTTTCCAACAATACGCTACTGCTGACGCGGTGCAAGCAGCTCGTATTCCCGAAGGGACCGACTTGGCCCAAGTCGCCCCCATCTTATGTGCTGGTATCACCGTCTACAAGGCTTTGAAGTCTGCTAACTTAATGGCCGGTCACTGGGTTGCTATCTCCGGTGCTGCTGGTGGTCTAGGTTCTTTAGCTGTTCAATACGCCAAGGCTATGGGTTACAGAGTCTTGGGTATTGACGGTGGTGAAGGTAAGGAAGAATTGTTCAGATCCATCGGTGGTGAAGTATTCATTGACTTCACTAAGGAAAAGGACATTGTCGGTGCTGTTCTAAAGGCAACTGACGGTGGTGCTCACGGTGTCATCAACGTTTCCGTTTCCGAAGCAGCTATTGAAGCTTCTACCAGATACGTTAGAGCTAACGGTACCACCGTTTTGGTAGGTATGCCAGCTGGTGCCAAGTGTTGTTCTGATGTCTTCAACCAAGTAGTCAAGTCCATCTCTATTGTTGGTTCTTACGTCGGTAACAGAGCTGACACAAGAGAAGCTTTGGACTTCTTCGCCAGAGGTTTGGTCAAGTCTCCAATCAAGGTTGTCGGATTGTCTACCTTGCCAGAAATTTACGAAAAGATGGAAAAGGGTCAAATCGTTGGTAGATACGTTGTTGACACATCTAAATAA")))

    val (markers, gaps) = Markers.getUniqueAggregatedKmersAndGaps(genes, Array(Markers.kmerCounts(genes, 21)), 21, false)
    val allgaps = gaps.flatten.toSet
    val reads = markers.zipWithIndex.map{ case (mG,i) => mG.map(m => Markers.genKmers(m.seq, 21).filter(!allgaps.contains(_)).map( k => Fasta.Entry("%s:%d:%d".format(genes(i).description, m.index, k.index), k.seq).toFastqEntry)).flatten}.flatten

    val quantifier = Quantification.Quantifier(markers.zipWithIndex.map{case (g, i) => g.map( k => Fasta.Entry("%s:%d".format(genes(i).description, k.index), k.seq))}.flatten, gaps.toArray.flatten.map(_.toFastaEntry), 21, false, false, minQual, Array.empty[String].toSet)

    reads.foreach{ r =>
      val res = quantifier.search(Array(r, r.revcomp))
      if (res.map(_.split(':')(0)).toSet.size > 1) {
        println("######\n>%s".format(r.description))
        res.map(println)
      }
    }

    genes.foreach{ f =>
      val res = quantifier.search(Array(f.toFastqEntry, f.revcomp.toFastqEntry))
      if (res.map(_.split(':')(0)).toSet.size > 1) {
        println("######\n>%s".format(f.description))
        res.map(println)
      }

    }
  }

}

}
