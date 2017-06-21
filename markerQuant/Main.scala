package markerQuant {

object Main{

  val ops = Map("markers" -> UniqueMarkers,
                //"quant"   -> Quantify,
                "test"    -> Test)

  def main(args: Array[String]) = {

    runAction(args)

  }

  def runAction(args: Array[String]) = {

    if (args.length < 1 || args(0).toLowerCase == "help" || !(ops contains args(0).toLowerCase)) {
      usage
    } else {

      val action = ops(args(0).toLowerCase)
      if ( args.length == 1  || args(1).toLowerCase == "help") {
        action.usage
      } else {
        action.main(args.slice(1, args.length))
       }
    }

  }

  def usage() = {

    println("markerQuant: Pipeline to identify markers")
    println("")
    println("Usage: panalysis <action> <action options>")
    println("Actions")
    ops.keys.toList.sorted.foreach{ k =>
      val action = ops(k)
      println("  %s%s%s -> %s".format(Console.BOLD, k, Console.RESET, action.description))
    }
  }


}


}
