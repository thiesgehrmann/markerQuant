import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult

import scala.collection.JavaConverters._

package markerQuant {

object Quantification {

  class Quantifier(val trees: Array[AhoCorasick], val markers: Array[Map[String,Markers.Kmer]], k: Int) {

    def searchTree(read: Fastq.Entry, treeID: Int) = {
      val results = this.trees(treeID).search(read.sequence.toString.getBytes).asInstanceOf[java.util.Iterator[SearchResult]].asScala.toArray.map(sr => sr.asInstanceOf[SearchResult].getOutputs.toArray.map(r => r.asInstanceOf[String])).flatten

      val aggmarkers = results.map{ mid =>  mid.substring(0,mid.lastIndexOf(':')) }.toSet
      val genes      = aggmarkers.map{ amid => amid.substring(0,amid.indexOf(':'))}.toSet

      val kmerCounts = Utils.CountMap(results.map( this.markers(treeID)(_).seq.toString))

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
        Utils.warning("ERROR! Read pair gave multiple targets!\n%s".format(f.map( r => "ID: %s\nNT: %s".format(r.description, r.sequence.toString)).mkString("\n")))
        Array("multitarget")
      } else if (targets.size == 0) {
        Array("unmapped")
      } else {
        markers
      }

    }
  }

  object Quantifier {

    def makeSingleTree(fastaKmers: Iterable[Fasta.Entry], fastaGaps: Iterable[Fasta.Entry], k: Int) = {

      val gaps = fastaGaps.zipWithIndex.map{ case (f,i) => Markers.Kmer(f.sequence, i) }.toSet

      val tree = new AhoCorasick()

      val markers = fastaKmers.map{f =>
        Markers.genKmers(f.sequence, k).zipWithIndex.filter{ case (s,i) =>
          //if (gaps.contains(s)) {
          //  Utils.message("Skipping %s becauase it is a gap!".format(s.seq.toString))
          //}
          !gaps.contains(s)
        }.map{ case (s,i) =>
          //Utils.message("Adding '%s:%d -> %s' to tree".format(f.description, i, s.seq.toString))
          "%s:%d".format(f.description, i) -> s
        }
      }.flatten.toMap

      markers.foreach{ case (k,v) =>
        tree.add(v.seq.toString.getBytes, k)
      }

      tree.prepare()
      (tree, markers)

    }

    def apply(fastaKmers: Iterable[Fasta.Entry], fastaGaps: Iterable[Fasta.Entry], k: Int, strandSpecific: Boolean, paired: Boolean) = {

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

      new Quantification.Quantifier(trees, markers, k)
    }

  }


}

}
