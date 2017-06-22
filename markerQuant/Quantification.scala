import org.arabidopsis.ahocorasick.AhoCorasick
import org.arabidopsis.ahocorasick.SearchResult

import scala.collection.JavaConverters._

package markerQuant {

object Quantification {

  def genTree(fastaKmers: Iterable[Fasta.Entry], fastaGaps: Iterable[Fasta.Entry], k: Int, strandSpecific: Boolean) = {
   
    val gaps = if (!strandSpecific) {
        fastaGaps.zipWithIndex.map{ case (f,i) => Markers.Kmer(f.sequence, i) }.toSet union fastaGaps.zipWithIndex.map{ case (f,i) => Markers.Kmer(f.sequence.revcomp, i) }.toSet
      } else {
        fastaGaps.zipWithIndex.map{ case (f,i) => Markers.Kmer(f.sequence, i) }.toSet
      }
 
    val tree = new AhoCorasick()

    fastaKmers.foreach{f =>
      Markers.genKmers(f.sequence, k).filter( s => !gaps.contains(s)).zipWithIndex.foreach{ case (s,i) =>
        tree.add(s.seq.toString.getBytes, "%d:%s".format(i, f.description))
      }
      if (!strandSpecific){
        Markers.genKmers(f.sequence.revcomp, k).filter( s => !gaps.contains(s)).zipWithIndex.foreach{ case (s,i) =>
          tree.add(s.seq.toString.getBytes, "RC_%d:%s".format(i, f.description))
        }
      }
    }

    tree.prepare()
    tree

  }

  def findInTree(tree: AhoCorasick, f: Array[Fastq.Entry]) = {
    /* Perform the search for each read end */
    val result = f.map(r => tree.search(r.sequence.toString.getBytes))
    /*For each result, remove the id at the start, and use only the stuff after the ':' */
    result.map( x => x.asInstanceOf[SearchResult].getOutputs.toArray.map(y => y.asInstanceOf[String])).flatten.map(m => m.substring(m.indexOf(':')+1)).toArray
  }

}

}
