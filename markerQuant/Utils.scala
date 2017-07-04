// Written by Thies Gehrmann thiesgehrmann@gmail.com 2015-2016

import scala.Console

package markerQuant {

///////////////////////////////////////////////////////////////////////////////

object Utils {

  def message(msg: String, messageType: String = "", ln: Boolean = true, or: Boolean = false) = { msg.split("\n").foreach(l => Console.err.print("%s%s%s".format(messageType, l, if (ln) "\n" else ""))) }
  def warning(msg: String, ln: Boolean =true) = message(msg, "WARNING: ", ln)
  def error(msg: String, ln: Boolean =true) = message(msg, "ERROR: ", ln)



  def suckLessZip[A](L: List[List[A]]): List[List[A]] = {
    /****************************************************************************
     * Zip together an arbitrary number of lists                                *
     ****************************************************************************/
  
    def zipHelper[A](previous: List[List[A]], current: List[A]): List[List[A]] = {
      val n_items = current.length
      var next = new collection.mutable.MutableList[List[A]]
      for ( i <- 0 until n_items) {
        next += previous(i) :+ current(i)
      }
      return next.toList
    }
    var previous  = L(0).map( x => List(x))
    for ( i <- 1 until L.length) {
      previous = zipHelper(previous, L(i));
    }
    return previous
  
  }

  def sumDicts[T](map1: Map[T,Int], map2: Map[T,Int]) : Map[T,Int] = {
    map1 ++ map2.map{ case (k,v) => k -> (v + map1.getOrElse(k,0)) }
  }

  class CountMap[T] {
    var counts = scala.collection.mutable.HashMap.empty[T,Int]

    def add(items: Iterable[T]){
       items.foreach{ i =>
         counts += (i -> (counts.getOrElse(i,0)+1))
       }
    }

    def get = this.counts
    def getOrElse(key: T, value: Int) = this.counts.getOrElse(key,value)
    def apply(key: T) = this.counts.getOrElse(key,0)
  }

  object CountMap {
    def apply[T](objects: Iterable[T]) = {
      val cm = new CountMap[T]
      cm.add(objects)
      cm
    }

    def empty[T] = new CountMap[T]
  }


}

}
