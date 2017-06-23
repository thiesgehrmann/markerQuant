import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult

import scala.collection.JavaConverters._

package markerQuant {

object Quantification {

  def genTree(fastaKmers: Iterable[Fasta.Entry], fastaGaps: Iterable[Fasta.Entry], k: Int) = {
  
    val gaps = fastaGaps.zipWithIndex.map{ case (f,i) => Markers.Kmer(f.sequence, i) }.toSet
 
    val tree = new AhoCorasick()

    fastaKmers.foreach{f =>
      Markers.genKmers(f.sequence, k).filter( s => !gaps.contains(s)).zipWithIndex.foreach{ case (s,i) =>
        //Utils.message("Adding '%s[%d] : %s' to tree".format(f.description, i, s.seq.toString))
        tree.add(s.seq.toString.getBytes, "%s".format(f.description))
      }
    }

    tree.prepare()
    tree

  }

  def findInTree(trees: Array[AhoCorasick], f: Array[Fastq.Entry]) = {
    /* Perform the search for each read end */
    //Utils.message("Checking sequences:\n%s\n".format(f.map(_.sequence.toString).mkString("\n")))
    val result = f.zipWithIndex.map{ case (r,i) => trees(i).search(r.sequence.toString.getBytes).asInstanceOf[java.util.Iterator[SearchResult]].asScala.toArray}.flatten
    //.zipWithIndex.map{ case (r, i) => // To test if the different pairs are matched with the tree. CHECK AGAIN FOR STRAND SPECIFIC!!!!
    //  println("%d: %d results".format(i, r.map(_.asInstanceOf[SearchResult].getOutputs.toArray).flatten.length))
    //  r
    //}.flatten
    result.map( x => x.asInstanceOf[SearchResult].getOutputs.toArray.map{y => 
      //Utils.message(y.asInstanceOf[String])
      y.asInstanceOf[String]
    }).flatten.toArray
  }

}

}
