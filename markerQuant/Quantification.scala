import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult

import scala.collection.JavaConverters._

package markerQuant {

object Quantification {

  class Quantifier(val trees: Array[AhoCorasick], val markers: Array[Map[String,Markers.Kmer]], k: Int, minQual: Byte) {

    def searchTree(read: Fastq.Entry, treeID: Int) = {
      val results = this.trees(treeID).search(read.sequence.toString.getBytes).asInstanceOf[java.util.Iterator[SearchResult]].asScala.toArray.map{ sr =>
        val result  = sr.asInstanceOf[SearchResult]
        val outputs = result.getOutputs.toArray.map(r => r.asInstanceOf[String])
        val lastIndex = result.getLastIndex
        (outputs, lastIndex)
      }

      val qualResults = results.filter{ case (outputs, lastIndex) =>
        !(read.quality.slice(lastIndex-k, lastIndex).exists(_ < this.minQual))
      }.map( _._1 ).flatten

      val aggmarkers = if ((qualResults.length == 0) && (results.length > 0)) {
        // We had hits, but they were of too low quality, add this to the count
        Array(Quantifier.lowQualRead).toSet
      } else {
        qualResults.map{ mid =>
        //  Take only the PROTID:MARKERLOC pair, and nothing else
          val idx = mid.indexOf(':', mid.indexOf(':')+1)
          if (idx >= 0) {
             mid.substring(0,idx)
          } else { mid }
        }.toSet
      }
      val genes = aggmarkers.map{ amid =>
        val idx = amid.lastIndexOf(':')
        if (idx >= 0) {
          amid.substring(0,amid.indexOf(':'))
        } else { amid }
      }.toSet

      if (genes.size > 1){
        Utils.message("MULTIMAPPED: \n%s\n%s\n".format(read.description, read.sequence.toString))
        qualResults.foreach( m => Utils.message("MULTIMAPPED: %s -> %s".format(m, this.markers(treeID)(m))))
        println("################")
      }

      if ((genes.size == 1) && genes.contains("augustus_masked-CENPK113-7D_chr_10_polished-processed-gene-4.26-mRNA-1")){
        Utils.message("Found TDH2: \n%s\n%s\n".format(read.description, read.sequence.toString))
        qualResults.foreach( m => Utils.message("FOUND TDH2: %s -> %s".format(m, this.markers(treeID)(m))))
        println("################")
      }

      //val kmerCounts = Utils.CountMap(results.map( this.markers(treeID)(_).seq.toString))

      //if (genes.size > 1){
      //  Utils.error("Read gave multiple targets:\n%s\nmarkers: %s".format(read.sequence.toString, results.map( m => "%s -> %s".format(m, this.markers(treeID)(m).seq.toString)).mkString("\n")) )
      //  kmerCounts.get.foreach{case (kmer, count) => println("%s: %d".format(kmer, count))}
      //}
      aggmarkers
    }

    def search(f: Array[Fastq.Entry]) : Array[String] = {
      /* Perform the search for each read end */
      //println("Checking pair:")
      val results = f.zipWithIndex.map{ case (r,i) => 
        this.searchTree(r, i)
      }
      //.zipWithIndex.map{ case (r, i) => // To test if the different pairs are matched with the tree. CHECK AGAIN FOR STRAND SPECIFIC!!!!
      //  println("%d: %d results".format(i, r.map(_.asInstanceOf[SearchResult].getOutputs.toArray).flatten.length))
      //  r
      //}.flatten
      val markers = results.flatten.toSet.toArray

      val targets = markers.map(_.split(':')(0)).toSet

      if (targets.size > 1){
        //Utils.warning("ERROR! Read pair gave multiple targets!\n%s".format(f.map( r => "ID: %s\nNT: %s".format(r.description, r.sequence.toString)).mkString("\n")))
        Array(Quantifier.multiTargetRead)
      } else if (targets.size == 0) {
        Array(Quantifier.unmappedRead)
      } else {
        markers
      }

    }
  }

  object Quantifier {

    val multiTargetRead = "__multitarget:0"
    val unmappedRead    = "__unmapped:0"
    val lowQualRead     = "__lowqual:0"
    val multiTargetReadTarget = "__multitarget"
    val unmappedReadTarget    = "__unmapped"
    val lowQualReadTarget     = "__lowqual"

    def makeSingleTree(fastaKmers: Iterable[Fasta.Entry], fastaGaps: Iterable[Fasta.Entry], k: Int) = {
      Utils.message("Making a tree")

      val gaps = fastaGaps.zipWithIndex.map{ case (f,i) => Markers.Kmer(f.sequence, i) }.toSet

      val tree = new AhoCorasick()

      val markers = fastaKmers.zipWithIndex.map{ case (f, fi) =>
        Utils.message("  \rAdding kmers for sequence %d ".format(fi), ln=false)
        Markers.genKmers(f.sequence, k).zipWithIndex.filter{ case (s,i) =>
          //if (gaps.contains(s)) {
          //  Utils.message("Skipping %s becauase it is a gap!".format(s.seq.toString))
          //}
          !gaps.contains(s)
        }.map{ case (s,i) =>
          tree.add(s.seq.toString.getBytes, "%s:%d".format(f.description, i))
          
//          Utils.message("Adding (%d) '%s:%d -> %s' to tree".format(fi, f.description, i, s.seq.toString))
          "%s:%d".format(f.description, i) -> s
        }
      }.flatten.toMap

      Utils.message("\rInserted %d sequences into tree\n".format(markers.size))

      tree.prepare()

      (tree, markers)

    }

    def apply(fastaKmers: Iterable[Fasta.Entry], fastaGaps: Iterable[Fasta.Entry], k: Int, strandSpecific: Boolean, paired: Boolean, minQual: Byte) = {
      Utils.message("Generating quantifier")

      def revcompHelper(sequences: Iterable[Fasta.Entry]) = { sequences.map( s => Fasta.Entry(s.description + ":RC", s.sequence.revcomp)) }

      val (tree, marker) = if (strandSpecific) { Quantifier.makeSingleTree(fastaKmers, fastaGaps, k) } else { Quantifier.makeSingleTree(fastaKmers ++ revcompHelper(fastaKmers), fastaGaps ++ revcompHelper(fastaGaps), k) }
      val (trees, markers) = if (!strandSpecific && paired) {
          (Array(tree, tree), Array(marker, marker))
        } else if (strandSpecific && paired) {
          val (secondTree, secondMarker) = Quantifier.makeSingleTree(revcompHelper(fastaKmers), revcompHelper(fastaGaps), k)
          (Array(tree, secondTree), Array(marker, secondMarker))
        } else if (!strandSpecific && paired) {
          (Array(tree, tree), Array(marker, marker))
        } else {
          (Array(tree, tree), Array(marker, marker))
        }

      val q = new Quantification.Quantifier(trees, markers, k, minQual)
      Utils.message("Added markers")

      q
    }

  }


}

}
